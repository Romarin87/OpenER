# OpenER

基于 ORCA 的化学基元反应过渡态工作流。该代码假定已准备好 TS 初始结构（xyz 或 inp），并且 ORCA/依赖库均已安装。

## 功能模块
- TS 几何优化：`! B3LYP D3BJ def2-SVP OptTS Freq`，遇到 SCF/步数问题自动重启并启用更紧的设置。
- 一阶鞍点验证：解析 ORCA 输出频率，要求虚频数量=1 且 |ν|>20 cm⁻¹。
- SOAP 去重：使用 DScribe 生成全局 SOAP 描述符，SQLite 持久化，cKDTree 搜索，阈值 S>0.999 或 d<1e-3 判定重复。
- IRC 路径：`! B3LYP D3BJ def2-SVP IRC`，Direction Both，MaxIter 50，StepSize 0.1，提取两端最低能几何作为 Reactant/Product guess。
- 端点最优化：与 TS 同级别 `! B3LYP D3BJ def2-SVP Opt Freq`，确认无虚频。
- Canonical SMILES 对比：OpenBabel 生成 SMILES，验证 IRC 端点与最优化结构拓扑一致。

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
- `TS_opt/`：ts_opt.inp/out 及 ts_opt_opt.xyz
- `IRC/`：irc.inp/out 及端点猜测 `irc_reactant_guess.xyz`, `irc_product_guess.xyz`
- `RP_opt/`：reactant.* 和 product.* 的优化/频率输出及 *_opt.xyz
- `run.log` 与 `result.json`：该结构的流程日志与摘要
可选 `--json summary.json` 仍可导出全局汇总。

## 关键文件
- `opener_workflow/config.py`：ORCA 关键字、SOAP 阈值、虚频判据。
- `opener_workflow/orca_runner.py`：ORCA 提交与自动重启，解析最终几何。
- `opener_workflow/analysis.py`：频率解析与鞍点/极小点判定。
- `opener_workflow/dedup.py`：SOAP 指纹计算、SQLite 存储与重复检测。
- `opener_workflow/irc.py`：IRC 输出解析，选择反应/生成物猜测结构。
- `opener_workflow/smiles_check.py`：Canonical SMILES 生成与匹配。
- `opener_workflow/pipeline.py`：整体流程（TS 优化 → 验证 → 去重 → IRC → 端点优化 → SMILES 校验）及 CLI。

## 说明
- 默认理论水平和网格可在 `config.py` 中调整。
- 输入 xyz 支持多个帧；inp 会自动解析 `* xyz` 块。
- 若 IRC 或优化失败，会在结果中标记 `failed`，详细见对应 `.out`。
