# OpenER

基于 ORCA 的化学基元反应过渡态工作流。该代码假定已准备好 TS 初始结构（xyz 或 inp），并且 ORCA/依赖库均已安装。

## 功能模块
- TS 几何优化：统一方法关键词（`method_keywords`），TS 用 `OptTS Freq`，`MaxIter` 默认沿用 ORCA（不写）；若重跑，仍用相同关键词，可选每隔 `ts_recalc_hess` 步重算 Hessian。
- 一阶鞍点验证：解析 ORCA 输出频率，要求虚频数量=1 且 |ν|>20 cm⁻¹。
- SOAP 去重：使用 DScribe 原子级 SOAP 描述符（`average="off"`），Laplacian 平均核，相似度阈值默认 0.99；指纹持久化 SQLite。若 IRC/SMILES 重试仍失败，会将该 TS 以状态写入 DB，下次命中直接跳过。
- IRC 路径：`IRC`，默认双向、`irc_maxiter=50`。若 SMILES 不一致，会按不匹配方向单向重跑，步数递增，重试次数由 `irc_max_retries` 控制；仍不符则失败并缓存。
- 端点最优化：与 TS 同级别 `Opt Freq`，确认无虚频。
- Canonical SMILES 对比：OpenBabel 生成 SMILES，验证 IRC 端点与最优化结构拓扑一致。
- 并行：仅在设置了 `OrcaSettings.launcher`（如 `srun`）时生效，可用 `PipelineConfig.max_workers` 控制并发处理多个 TS 输入。

## 快速开始
```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt

# 运行工作流
python -m opener_workflow.pipeline \
  --ts-dir path/to/ts_xyz_or_inp \
  --db data/soap_db.sqlite \
  --charge 0 --mult 1
```

未指定 `--workdir` 时默认生成 `runs/<时间戳>/`，避免覆盖旧结果；如需固定路径可显式传入。运行结束后，每个结构的结果放在 `workdir/<结构标签>/`，并包含：
- `TS_opt/`：ts_opt.inp/out 及 ORCA 写出的 ts_opt.xyz（最终结构）。
- `IRC/`：irc.inp/out 及 ORCA 写出的端点 `irc_IRC_B.xyz`（Backward）和 `irc_IRC_F.xyz`（Forward）。
- `RP_opt/`：reactant.* 和 product.* 的 Opt+Freq 输出及 ORCA 写出的 reactant.xyz / product.xyz。
- `run.log` 与 `result.json`：该结构的流程日志与摘要
可选 `--json summary.json` 仍可导出全局汇总。

## 关键文件
- `opener_workflow/config.py`：ORCA 关键字、SOAP 阈值、虚频判据。
- `opener_workflow/orca_runner.py`：ORCA 提交与自动重启，直接读取 ORCA 输出的 xyz（包括 TS、IRC 端点、端点优化）。
- `opener_workflow/analysis.py`：频率解析、鞍点/极小点判定，以及 SMILES 生成/对比。
- `opener_workflow/dedup.py`：SOAP 指纹计算、SQLite 存储与重复检测。
- `opener_workflow/pipeline.py`：整体流程（TS 优化 → 验证 → 去重 → IRC → 端点优化 → SMILES 校验）及 CLI。

## 说明
- 默认理论水平和网格可在 `config.py` 中调整。
- 输入 xyz 支持多个帧；inp 会自动解析 `* xyz` 块。
- 若 IRC 或优化失败，会在结果中标记 `failed`，详细见对应 `.out`。
