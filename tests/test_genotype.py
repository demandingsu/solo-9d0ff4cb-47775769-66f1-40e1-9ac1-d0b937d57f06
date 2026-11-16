"""基因型（GT）扩展测试：单样本上传、GT 归一、冲突失败、版本切换、旧库迁移。"""
import json
import sqlite3
from pathlib import Path

import pytest
from fastapi.testclient import TestClient

SAMPLES = Path(__file__).resolve().parent.parent / "samples"

# 升级前（无 GT 功能）的旧库结构，用于验证旧库存量迁移与“未知”标注
LEGACY_SCHEMA = """
CREATE TABLE datasets (
    id               INTEGER PRIMARY KEY AUTOINCREMENT,
    name             TEXT NOT NULL,
    version          INTEGER NOT NULL DEFAULT 0,
    published_job_id INTEGER,
    created_at       TEXT NOT NULL,
    updated_at       TEXT NOT NULL
);
CREATE TABLE jobs (
    id            INTEGER PRIMARY KEY AUTOINCREMENT,
    dataset_id    INTEGER NOT NULL REFERENCES datasets(id),
    status        TEXT NOT NULL CHECK (status IN ('ready', 'failed')),
    fasta_name    TEXT NOT NULL,
    vcf_name      TEXT NOT NULL,
    chrom         TEXT,
    total_records INTEGER NOT NULL,
    bad_count     INTEGER NOT NULL,
    errors_json   TEXT NOT NULL DEFAULT '[]',
    created_at    TEXT NOT NULL
);
CREATE TABLE job_variants (
    job_id     INTEGER NOT NULL REFERENCES jobs(id),
    chrom      TEXT NOT NULL,
    pos        INTEGER NOT NULL,
    ref        TEXT NOT NULL,
    alt        TEXT NOT NULL,
    lines_json TEXT NOT NULL
);
CREATE TABLE variants (
    dataset_id INTEGER NOT NULL REFERENCES datasets(id),
    chrom      TEXT NOT NULL,
    pos        INTEGER NOT NULL,
    ref        TEXT NOT NULL,
    alt        TEXT NOT NULL,
    lines_json TEXT NOT NULL
);
"""


@pytest.fixture()
def client(tmp_path, monkeypatch):
    monkeypatch.setenv("DATABASE_PATH", str(tmp_path / "test.db"))
    from app.main import app

    with TestClient(app) as c:
        yield c


def create_dataset(client, name="ds") -> int:
    r = client.post("/datasets", json={"name": name})
    assert r.status_code == 201, r.text
    return r.json()["id"]


def upload(client, dataset_id: int, vcf_name: str):
    with open(SAMPLES / "chr1.fa", "rb") as fa, open(SAMPLES / vcf_name, "rb") as vcf:
        return client.post(
            f"/datasets/{dataset_id}/uploads",
            files={"fasta": ("chr1.fa", fa, "text/plain"), "vcf": (vcf_name, vcf, "text/plain")},
        )


def upload_text(client, dataset_id: int, vcf_text: str):
    return client.post(
        f"/datasets/{dataset_id}/uploads",
        files={
            "fasta": ("chr1.fa", (SAMPLES / "chr1.fa").read_bytes(), "text/plain"),
            "vcf": ("inline.vcf", vcf_text.encode(), "text/plain"),
        },
    )


def publish(client, dataset_id: int, job_id: int, expected_version: int):
    return client.post(
        f"/datasets/{dataset_id}/publish",
        json={"job_id": job_id, "expected_version": expected_version},
    )


def gt_by_locus(variants):
    return {(v["pos"], v["ref"], v["alt"]): v["gt"] for v in variants}


def test_single_sample_upload_and_publish(client):
    """单样本 VCF：GT 随位点归并，1/0 归一为 0/1，./. 标为未知，发布后随版本可查。"""
    did = create_dataset(client)

    r = upload(client, did, "sample1.vcf")
    assert r.status_code == 201, r.text
    job = r.json()
    assert job["status"] == "ready"
    assert job["sample_name"] == "NA001"
    assert job["total_records"] == 6 and job["bad_count"] == 0

    gts = gt_by_locus(job["variants"])
    assert gts == {
        (4, "G", "GTA"): "0/1",   # 行 6 (0/1) 与行 7 (1/0) 归并，1/0 归一为 0/1
        (4, "GTA", "G"): "1/1",
        (17, "A", "G"): "0/1",
        (18, "T", "G"): None,     # 源行 ./. -> 未知，不推定为 0/0
        (20, "G", "GA"): "0/0",   # 显式 0/0 原样保留
    }

    # 发布前查询仍是旧版（无样本、无位点）
    detail = client.get(f"/datasets/{did}").json()
    assert detail["version"] == 0 and detail["sample_name"] is None

    r = publish(client, did, job["id"], expected_version=0)
    assert r.status_code == 200, r.text

    # 发布后：样本名、GT、位点同属 v1
    detail = client.get(f"/datasets/{did}").json()
    assert detail["version"] == 1
    assert detail["sample_name"] == "NA001"
    assert gt_by_locus(detail["variants"]) == gts
    # 列表接口同样带样本名
    listed = client.get("/datasets").json()
    assert listed[0]["sample_name"] == "NA001"


def test_no_sample_vcf_gt_unknown(client):
    """无样本 VCF 继续接受：所有位点 GT 明确为未知（null），不推定为 0/0。"""
    did = create_dataset(client)
    job = upload(client, did, "batch1.vcf").json()
    assert job["status"] == "ready" and job["sample_name"] is None
    assert all(v["gt"] is None for v in job["variants"])

    assert publish(client, did, job["id"], 0).status_code == 200
    detail = client.get(f"/datasets/{did}").json()
    assert detail["sample_name"] is None
    assert detail["variant_count"] == 5
    assert all(v["gt"] is None for v in detail["variants"])


def test_gt_conflict_fails_job_and_keeps_old_version(client):
    """归并到同一位点的 GT 冲突 -> 整项任务失败，不可发布，旧版不变。"""
    did = create_dataset(client)
    job1 = upload(client, did, "sample1.vcf").json()
    assert publish(client, did, job1["id"], 0).status_code == 200
    before = client.get(f"/datasets/{did}").json()

    r = upload(client, did, "gt_conflict.vcf")
    assert r.status_code == 201, r.text
    job2 = r.json()
    assert job2["status"] == "failed"
    assert job2["total_records"] == 3 and job2["bad_count"] == 1
    err = job2["errors"][0]
    assert err["line"] == 7  # 后读行（物理行 7）被判冲突，不会覆盖先读行
    assert "GT 冲突" in err["reason"] and "0/1" in err["reason"] and "1/1" in err["reason"]

    # 失败任务发布 -> 422，已发布数据不变
    r = publish(client, did, job2["id"], expected_version=1)
    assert r.status_code == 422
    after = client.get(f"/datasets/{did}").json()
    assert after["version"] == before["version"] == 1
    assert after["sample_name"] == before["sample_name"] == "NA001"
    assert gt_by_locus(after["variants"]) == gt_by_locus(before["variants"])


def test_publish_switches_sample_and_gt_with_version(client):
    """发布时样本名、GT 与位点随版本一起切换；查询只返回同一已提交版本的数据。"""
    did = create_dataset(client)
    job1 = upload(client, did, "sample1.vcf").json()
    assert publish(client, did, job1["id"], 0).status_code == 200

    job2 = upload(client, did, "sample2.vcf").json()
    assert job2["sample_name"] == "NA002"
    assert publish(client, did, job2["id"], 1).status_code == 200

    detail = client.get(f"/datasets/{did}").json()
    assert detail["version"] == 2
    assert detail["sample_name"] == "NA002"
    gts = gt_by_locus(detail["variants"])
    assert gts[(4, "G", "GTA")] == "1/1"   # 来自 NA002，不再是 NA001 的 0/1
    assert gts[(17, "A", "G")] == "1/1"
    assert gts[(18, "T", "G")] == "0/1"

    # 再发布无样本 VCF：样本名与 GT 一起切回“未知”
    job3 = upload(client, did, "batch1.vcf").json()
    assert publish(client, did, job3["id"], 2).status_code == 200
    detail = client.get(f"/datasets/{did}").json()
    assert detail["version"] == 3
    assert detail["sample_name"] is None
    assert all(v["gt"] is None for v in detail["variants"])


def test_invalid_gt_and_missing_gt_are_bad_rows(client):
    """非法 GT、FORMAT 缺 GT、样本列缺值 -> 坏行，任务失败。"""
    did = create_dataset(client)
    vcf = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA001\n"
        "chr1\t17\t.\tA\tG\t50\tPASS\t.\tGT\t0|1\n"      # 行 3： phased 不在支持集
        "chr1\t17\t.\tA\tG\t50\tPASS\t.\tGT\t1/2\n"      # 行 4： 非二倍体支持集
        "chr1\t17\t.\tA\tG\t50\tPASS\t.\tDP\t30\n"       # 行 5： FORMAT 无 GT
        "chr1\t17\t.\tA\tG\t50\tPASS\t.\tDP:GT\t30\n"    # 行 6： 样本列值数不足，缺 GT 值
        "chr1\t20\t.\tG\tGA\t50\tPASS\t.\tGT:DP\t1/1:30\n"  # 行 7： 多子字段正常读取
    )
    job = upload_text(client, did, vcf).json()
    assert job["status"] == "failed"
    assert job["total_records"] == 5 and job["bad_count"] == 4
    reasons = {e["line"]: e["reason"] for e in job["errors"]}
    assert "GT 非法" in reasons[3] and "GT 非法" in reasons[4]
    assert "FORMAT 缺少 GT" in reasons[5]
    assert "缺少 GT 值" in reasons[6]
    assert publish(client, did, job["id"], 0).status_code == 422


def test_multi_sample_vcf_rejected(client):
    """多样本 VCF 属文件级错误 -> 400，不产生任务。"""
    did = create_dataset(client)
    vcf = (
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tNA001\tNA002\n"
        "chr1\t17\t.\tA\tG\t50\tPASS\t.\tGT\t0/1\t1/1\n"
    )
    r = upload_text(client, did, vcf)
    assert r.status_code == 400
    assert "单样本" in r.json()["detail"]


def test_legacy_database_migrates_and_marks_gt_unknown(tmp_path, monkeypatch):
    """旧库存量（无 GT 列的库）：迁移后可查，GT 明确为未知，且可继续发布带 GT 的新版本。"""
    db = tmp_path / "legacy.db"
    conn = sqlite3.connect(db)
    conn.executescript(LEGACY_SCHEMA)
    conn.execute(
        "INSERT INTO datasets (name, version, published_job_id, created_at, updated_at) "
        "VALUES ('legacy', 1, 1, '2025-01-01T00:00:00+00:00', '2025-01-01T00:00:00+00:00')"
    )
    conn.execute(
        "INSERT INTO jobs (id, dataset_id, status, fasta_name, vcf_name, chrom, "
        "total_records, bad_count, errors_json, created_at) "
        "VALUES (1, 1, 'ready', 'ref.fa', 'old.vcf', 'chr1', 1, 0, '[]', '2025-01-01T00:00:00+00:00')"
    )
    conn.execute(
        "INSERT INTO variants (dataset_id, chrom, pos, ref, alt, lines_json) "
        "VALUES (1, 'chr1', 17, 'A', 'G', '[5]')"
    )
    conn.commit()
    conn.close()

    monkeypatch.setenv("DATABASE_PATH", str(db))
    from app.main import app

    with TestClient(app) as client:
        # 旧版本数据：GT 明确标为未知（null），不会推成 0/0
        detail = client.get("/datasets/1").json()
        assert detail["version"] == 1 and detail["sample_name"] is None
        assert detail["variants"][0]["gt"] is None
        assert detail["variants"][0]["lines"] == [5]

        # 旧库上继续完整流程：上传单样本 VCF 并按版本发布
        job = upload(client, 1, "sample1.vcf").json()
        assert job["status"] == "ready"
        r = publish(client, 1, job["id"], expected_version=1)
        assert r.status_code == 200, r.text
        detail = client.get("/datasets/1").json()
        assert detail["version"] == 2 and detail["sample_name"] == "NA001"
        assert gt_by_locus(detail["variants"])[(4, "G", "GTA")] == "0/1"
