"""上传内容处理：逐条校验、归一化、归并等价位点。"""
from __future__ import annotations

from dataclasses import dataclass

from .normalize import NormalizeError, normalize
from .parse import VcfRow, parse_fasta, parse_vcf_rows

# VCF REF/ALT 允许的碱基字符
VCF_BASES = frozenset("ACGTN")

# 可接受的二倍体 GT 原始写法（1/0 归一为 0/1；./. 表示未知）
GT_KNOWN = frozenset({"0/0", "0/1", "1/1"})
GT_MISSING = "./."


@dataclass
class RowError:
    """一条坏行的诊断信息。"""

    line: int
    reason: str
    chrom: str | None = None
    pos: str | None = None
    ref: str | None = None
    alt: str | None = None

    def to_dict(self) -> dict:
        return {
            "line": self.line,
            "reason": self.reason,
            "chrom": self.chrom,
            "pos": self.pos,
            "ref": self.ref,
            "alt": self.alt,
        }


@dataclass
class NormVariant:
    """归一化后的位点，lines 为归并到该位点的所有原始行号。

    gt 为该位点的二倍体基因型（0/1 与 1/0 已归一为 0/1）；
    None 表示未知（无样本列上传或源行 ./.），绝不推定为 0/0。
    """

    pos: int
    ref: str
    alt: str
    lines: list[int]
    gt: str | None = None


def _check_row(fasta_chrom: str, row: VcfRow) -> tuple[int, str, str]:
    """单行结构校验，返回 (pos, ref, alt)；不合法抛 NormalizeError。"""
    if row.parse_error:
        raise NormalizeError(row.parse_error)
    if row.chrom != fasta_chrom:
        raise NormalizeError(f"CHROM={row.chrom} 与 FASTA 染色体 {fasta_chrom} 不一致")
    try:
        pos = int(row.pos)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        raise NormalizeError(f"POS 不是整数：{row.pos!r}")
    ref = (row.ref or "").upper()
    alt = (row.alt or "").upper()
    if not ref or any(b not in VCF_BASES for b in ref):
        raise NormalizeError(f"REF 含非法字符：{row.ref!r}（仅支持 ACGTN）")
    if "," in alt:
        raise NormalizeError("ALT 含逗号：仅支持双等位变异")
    if not alt or any(b not in VCF_BASES for b in alt):
        raise NormalizeError(f"ALT 非法：{row.alt!r}（仅支持 ACGTN，不支持 <符号>/断点/*/.）")
    return pos, ref, alt


def _extract_gt(row: VcfRow) -> str | None:
    """从 FORMAT/样本列读取二倍体 GT。

    返回 "0/0" / "0/1" / "1/1"（1/0 已归一为 0/1）；
    ./. 或无样本列时返回 None（未知）。非法 GT 抛 NormalizeError（坏行）。
    """
    if row.format_field is None:
        return None  # 无样本列的上传，GT 一律未知
    fmt_keys = row.format_field.split(":")
    if "GT" not in fmt_keys:
        raise NormalizeError(f"FORMAT 缺少 GT 字段：{row.format_field!r}")
    idx = fmt_keys.index("GT")
    sample_parts = (row.sample_field or "").split(":")
    raw = sample_parts[idx].strip() if idx < len(sample_parts) else ""
    if not raw:
        raise NormalizeError("样本列缺少 GT 值")
    if raw == "1/0":
        return "0/1"  # 0/1 与 1/0 视为相同，归一为 0/1
    if raw in GT_KNOWN:
        return raw
    if raw == GT_MISSING:
        return None  # 显式缺失 -> 未知
    raise NormalizeError(f"GT 非法：{raw!r}（仅支持二倍体 0/0、0/1、1/0、1/1、./.）")


def process_upload(
    fasta_text: str, vcf_text: str
) -> tuple[str, str | None, int, list[NormVariant], list[RowError]]:
    """处理一次上传。

    返回:
        (chrom, sample_name, total_records, variants, errors)
        sample_name 为单样本 VCF 的样本列名，无样本 VCF 为 None；
        variants 已按 (pos, ref, alt) 归并等价位点并排序，lines 为升序原始行号，
        gt 为归并后的基因型（None=未知）。

    归并到同一位点的行若 GT 冲突（两个不同的已知 GT），冲突行记为坏行、
    整项任务失败，不会用后读行覆盖先读行；未知（./.）让位于已知，不算冲突。

    异常:
        ParseError: 文件级致命错误（由调用方转为 HTTP 400）。
    """
    chrom, seq = parse_fasta(fasta_text)
    rows, sample_name = parse_vcf_rows(vcf_text)

    merged: dict[tuple[int, str, str], list[int]] = {}
    gt_by_key: dict[tuple[int, str, str], str | None] = {}
    gt_line: dict[tuple[int, str, str], int] = {}  # 确立各位点 GT 的原始行号
    errors: list[RowError] = []
    for row in rows:
        try:
            pos, ref, alt = _check_row(chrom, row)
            npos, nref, nalt = normalize(seq, pos, ref, alt)
            gt = _extract_gt(row)
        except NormalizeError as exc:
            errors.append(
                RowError(line=row.line, reason=str(exc), chrom=row.chrom, pos=row.pos, ref=row.ref, alt=row.alt)
            )
            continue
        key = (npos, nref, nalt)
        merged.setdefault(key, []).append(row.line)
        if key not in gt_by_key:
            gt_by_key[key] = gt
            gt_line[key] = row.line
        elif gt is not None:
            if gt_by_key[key] is None:
                gt_by_key[key] = gt  # 未知（./.）让位于已知
                gt_line[key] = row.line
            elif gt != gt_by_key[key]:
                errors.append(
                    RowError(
                        line=row.line,
                        reason=(
                            f"GT 冲突：与第 {gt_line[key]} 行归并到同一位点 "
                            f"{npos} {nref}>{nalt}，GT {gt_by_key[key]} 与 {gt} 不一致"
                            f"（整项任务失败，不会用后读行覆盖先读行）"
                        ),
                        chrom=row.chrom,
                        pos=row.pos,
                        ref=row.ref,
                        alt=row.alt,
                    )
                )

    variants = [
        NormVariant(pos=p, ref=r, alt=a, lines=sorted(lines), gt=gt_by_key[(p, r, a)])
        for (p, r, a), lines in sorted(merged.items())
    ]
    return chrom, sample_name, len(rows), variants, errors
