# VCF 变异归一化服务

面向测序中心数据员的跨批次变异归一化 API：不同批次 VCF 对同一变异常有不同写法
（锚点位置、共有碱基、重复区内的坐标漂移），本服务将其归一化为唯一表示并归并等价位点，
供跨批次比对；单样本 VCF 同时读取基因型（GT），供下游按规范化位点查询。
FastAPI + SQLite，Docker Compose 一键启动。

## 功能与规则

**输入**：单染色体 FASTA + 仅含双等位 SNV / 插入 / 缺失的 VCF（无样本或单样本）。

**归一化流水线**（`app/normalize.py`，与 vt normalize / bcftools norm 左对齐语义一致）：

1. **逐条核对 REF**：VCF 每行的 REF 必须与上传 FASTA 对应位置完全一致（大小写不敏感）；
2. **消去共有碱基**：先修剪 REF/ALT 共有后缀，再修剪共有前缀；
3. **保留合法锚点**：纯插入/缺失保留一个锚点碱基（优先左锚点，染色体左端用右锚点）；
4. **左对齐**：串联重复中的插入/缺失左移至最左等价表示；
5. **归并等价位点**：归一化后 (POS, REF, ALT) 相同的行合并，保留**各原始行号**
   （VCF 文件物理行号，从 1 起，含头部行）。

**基因型（GT）**：

- 无样本 VCF（头行 8 列）继续接受，所有位点 GT 为 `null`（未知）；
- 单样本 VCF（头行 10 列：FORMAT + 一个样本列）读取 FORMAT 中的二倍体 GT，
  仅接受 `0/0`、`0/1`、`1/0`、`1/1`、`./.`，其中 **`1/0` 归一为 `0/1`**（两者视为相同）；
- `./.` 表示未知，输出为 `gt: null`；**未知绝不推定为 `0/0`**；
- 不同原始行归并到同一位点后，若已知 GT 冲突（如 `0/1` vs `1/1`），**整项暂存任务失败**，
  不会用后读行覆盖先读行；未知（`./.`）与已知归并时取已知，不算冲突；
- 多样本 VCF（头行 >10 列）为文件级错误（HTTP 400）。

**坏行**（任一即判坏行，含坏行的任务不可发布）：REF 与参考不符、CHROM 与 FASTA 不一致、
POS 非整数/越界、REF 超出染色体末端、多等位（ALT 含逗号）、符号/星号 ALT、MNP、
替换伴随长度变化的复杂变异、REF==ALT、列数不足、GT 非法/缺失、归并后 GT 冲突。

**文件级致命错误**（HTTP 400，不产生任务）：FASTA 含多条序列、非法碱基字符、
VCF 缺 `#CHROM` 头、多样本 VCF、非 UTF-8 文件。

## 一致性模型

- 上传生成**暂存任务**（job），处理期间查询仍返回**完整旧版**；
- 发布按 `数据集号 + expected_version`（查询所得版本）在**单事务内整体替换**，
  **样本名、GT 与位点随版本一起切换**，查询只返回同一已提交版本的数据；
- **并发替换只允许版本匹配者成功**（乐观锁）：不匹配返回 409 及当前版本；
- **坏行任务或发布失败不改动旧版**：失败任务发布返回 422 及逐行原因，任何失败均回滚；
- 旧版本库（无 GT 列）启动时自动迁移（`ALTER TABLE` 补列），**旧库存量的 GT 明确标为
  未知（`gt: null`），不会推成 0/0**。

## 快速开始

```bash
docker compose up --build
```

- API 文档（Swagger UI）：http://localhost:8000/docs
- 健康检查：http://localhost:8000/healthz
- SQLite 数据持久化在 named volume `variant-data`（容器内 `/data/variants.db`）

### 一键演示（等价变异归并 + 版本乐观锁 + 坏行隔离 + 单样本基因型）

```bash
python3 scripts/demo.py            # 默认 http://localhost:8000
```

### curl 示例

```bash
# 1. 创建数据集（version 从 0 起）
curl -s -X POST localhost:8000/datasets -H 'Content-Type: application/json' \
     -d '{"name":"my-dataset"}'

# 2. 上传 FASTA + VCF，生成暂存任务（不影响已发布数据）
curl -s -X POST localhost:8000/datasets/1/uploads \
     -F 'fasta=@samples/chr1.fa' -F 'vcf=@samples/batch1.vcf'

# 3. 查询当前版本（处理期间/发布前仍返回旧版）
curl -s localhost:8000/datasets/1

# 4. 按查询所得版本发布（乐观锁；不匹配则 409，坏行任务则 422）
curl -s -X POST localhost:8000/datasets/1/publish \
     -H 'Content-Type: application/json' -d '{"job_id":1,"expected_version":0}'

# 5. 查询任务状态（坏行明细 / 暂存结果）
curl -s localhost:8000/datasets/1/jobs/1
```

## API 一览

| 方法 | 路径 | 说明 |
|---|---|---|
| POST | `/datasets` | 创建数据集，返回 `version=0` |
| GET | `/datasets` | 列出全部数据集及版本（含当前样本名） |
| GET | `/datasets/{id}` | 查询当前已发布版本（坐标、REF/ALT、GT、各原始行号、样本名） |
| POST | `/datasets/{id}/uploads` | 上传 FASTA+VCF 并归一化，生成暂存任务（无样本或单样本 VCF） |
| GET | `/datasets/{id}/jobs/{job_id}` | 任务状态、坏行明细、暂存的归一化结果（含 GT） |
| POST | `/datasets/{id}/publish` | 按 `expected_version` 乐观锁发布；409=版本冲突，422=任务含坏行 |
| GET | `/healthz` | 健康检查 |

**GT 输出约定**：`gt` 为 `"0/0"` / `"0/1"` / `"1/1"`（`1/0` 已归一为 `"0/1"`）或 `null`；
`null` 明确表示**未知**（无样本列上传、源行 `./.`、旧版本数据），不会推定为 `0/0`。
`sample_name` 为当前已发布版本的样本列名，随版本一起切换，无样本时为 `null`。

## 演示样本（`samples/`）

`chr1.fa`：40 bp 演示染色体，5–12 位为 TA 串联重复（`GCAG TATATATA GCCG...`）。

`batch1.vcf` 与 `batch2.vcf` 是**同一批变异在两批数据中的不同写法**：

| 变异 | 批次 1 写法（行号） | 批次 2 写法（行号） | 归一化结果 |
|---|---|---|---|
| TA 重复区插入 TA | `6 A>ATA`(5)、`8 A>ATA`(6) | `10 A>ATA`(5)、`12 A>ATA`(6)、`5 T>TAT`(7) | **4 G>GTA** |
| TA 重复区缺失 TA | `6 ATA>A`(7) | `4 GTA>G`(8)、`5 TAT>T`(9) | **4 GTA>G** |
| SNV | `17 A>G`(8) | `17 A>G`(10) | **17 A>G** |
| SNV（带共有碱基写法） | `18 TT>GT`(9) | `18 T>G`(11) | **18 T>G** |
| 非重复区插入 | `20 G>GA`(10) | `20 G>GA`(12) | **20 G>GA** |

两批发布后归一化结果**完全一致**（仅原始行号不同），可直接跨批次比对。

`bad.vcf`：5 条坏行（REF 不符、多等位、MNP、越界、CHROM 不符），用于演示坏行隔离。

`sample1.vcf` / `sample2.vcf`：单样本 VCF（样本列 `NA001` / `NA002`），位点与批次写法对应，
演示 GT 读取（含 `1/0`→`0/1` 归一、`./.`→未知、显式 `0/0`）与随版本切换。

`gt_conflict.vcf`：同一位点（TA 重复区插入）两行 GT 分别为 `0/1` 与 `1/1`，
归并后冲突 → 整项任务失败，用于演示 GT 冲突隔离。

## 本地开发与测试

```bash
python3 -m venv .venv && .venv/bin/pip install -r requirements-dev.txt
.venv/bin/python -m pytest tests/ -q          # 34 个测试：算法 + API + 并发发布 + 基因型
DATABASE_PATH=data/variants.db .venv/bin/python -m uvicorn app.main:app --reload
```

## 目录结构

```
app/
  main.py        FastAPI 入口与路由
  normalize.py   归一化算法（核对/修剪/锚点/左对齐）
  parse.py       FASTA/VCF 结构解析（含单样本列识别）
  service.py     逐行校验、GT 读取与冲突检测、归并等价位点
  store.py       SQLite 存储、旧库迁移与乐观锁发布事务
  schemas.py     请求/响应模型
samples/         演示用 FASTA 与六份 VCF
scripts/demo.py  端到端演示脚本（仅标准库）
tests/           pytest 测试
Dockerfile / docker-compose.yml
```
