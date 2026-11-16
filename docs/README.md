# 文档索引

七份 LaTeX 文档，一份一个文件夹，`.tex` 与编译好的 `.pdf` 放在一起。
编译命令（`xelatex`，中文文档必须用它，`pdflatex` 会掉字）：

```
cd docs/ref_cn && latexmk -xelatex ref_uncertainty_cn.tex
```

| 文件夹 | 文档 | 语言 | 页数 | 内容 |
|---|---|---|---|---|
| `ref_cn/` | `ref_uncertainty_cn.tex` | 中文 | 10 | 不确定度与误差传递参考：测量链、C_p / C_D / θ / Re 四条传递链、呈现规则、代码映射、口试问答 |
| `ref_en/` | `ref_uncertainty_en.tex` | 英文 | 10 | 同上，逐节镜像，可直接引用进英文报告 |
| `ops_en/` | `lab_operation_sheet.tex` | 英文 | 6 | 现场操作流程：出发前清单、开机、±30 对中判据、Pitot 判定、run 矩阵、读数流程、离场四查、**记录项→CSV 列→cfg 参数对照**、故障表、时间预算、打勾清单 |
| `ops_cn/` | `lab_operation_sheet_cn.tex` | 中文 | 6 | 同上，中文版 |
| `prelab/` | `prelab_circular_cylinder.tex` | 英文 | 8 | 口试前预习：Cp 的测压管形式推导、势流与 D'Alembert、分离与平台、trip wire、Pitot 判别与 ±30 对齐、run plan、减速六步（含 Roshko 三式）、不确定度计划、bench 自检、**口试 11 问逐条答案**、logbook 清单 |
| `prelab/` | `prelab_onepager.tex` | 英文 | 1 | 上面那份的 logbook 压缩版：目标、器材数字、11 步流程、减速顺序、不确定度计划、要记的列、run 矩阵、11 问提示，带一个画孔位环的空框，可直接抄进 Immutalog |
| `during_lab/` | `during_lab_sheet.tex` | 英文 | 7 | 现场记录模板：session 头、室温三读、Pitot 识别、每 run 一张对齐表、首轮 5 项 checklist、每个 run 一页 30 行读数表（三列重复读数）、收尾核对、七个数抄录表 |
| `during_lab/` | `logbook_entries.tex` | 英文 | 4 | logbook 条目模板：每条 entry 的标题与字段（session 头、环境三读、Pitot 识别、±30 对齐判定）；每 run 一张卡片含 **28 孔 ×(r1,r2,r3) 读数表**（两栏 14 行、列名对齐 `readings.csv`）+ Did / Saw / Decided 三行；收尾一致性检查、汇总表、异常处理速查 |

## 编译注意

- 中文两份用 `ctexart`，必须 `xelatex`；英文两份是 `article` + `fontspec`，`xelatex` 或 `pdflatex` 都行。
- 图不内嵌，需要时从 `../results/` 取 `DEMO_*.pdf`（矢量），加 `\includegraphics[width=\linewidth]{}`。
- 编译产生的 `.aux/.log/.xdv/.fls/.fdb_latexmk/.out` 不应留在文件夹里，可随时删。

## 与代码的对应

参考文档与 prelab 里的符号、阈值、列名都和 `../cylinder_lab.m`、`../cpUncertainty.m`、`../runUncertainty.m`
里的一致；改代码里的阈值（`cfg.uRead`、`cfg.uTheta`、`cfg.alignWarn`、`cfg.sepUncK`）时，
`ref_cn` / `ref_en` 第 2、5、7 节的表格要同步改，否则文档和脚本会脱节。
