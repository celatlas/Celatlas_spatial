# Celatlas Spatial v1.8.0 发布说明

## 定位

v1.8.0 是相对 v1.7.0 的独立重大更新。建议在 GitHub 上使用单独的 `v1.8.0` 分支和 `v1.8.0` tag，并把 Linux x86_64 部署包上传到对应 Release。

## 主要变化

- 新增统一公共 CLI：`celatlas count`、`celatlas reanalyze`、`celatlas mkreport`。
- 新增 ST、SX、SN 工作流和 Python runner；支持配置/CSV 批处理。
- 新增 FASTQ 发现与输入预检、运行计划、失败报告、步骤级执行和断点续跑。
- 新增已有 count 结果的 reanalysis 与仅报告生成模式。
- 新增 StarDist 细胞分割入口及 cell-level 分析：细胞级 10X matrix、AnnData、空间坐标、QC、聚类和 marker 输出。
- 扩展 HE/ssDNA/gene expression 图像流程、报告模板和可选模型交付。
- 新增 Linux 发布包：Conda 环境定义、安装/验证脚本、双语文档和合成 FASTQ 冒烟数据。
- 保留 `Celatlas.sh`、`Celatlas_FFPE.sh`、`Celatlas_SN.sh` 和 `Celatlas_reanalysis.sh` 兼容入口。

## 版本与环境

- Python：3.11+
- Conda 环境：`celatlas18`
- 支持平台：Linux x86_64
- 细胞分割模型：发布包内提供可选 StarDist H&E 模型；Swin 模型按 `release/MISSING_RUNTIME_ASSETS.md` 准备。

## 发布资产

- `celatlas-spatial-v1.8.0-linux-x86_64.tar.gz`：完整 Linux 部署包
- `celatlas_spatial-1.8.0-py3-none-any.whl`：Python wheel（包含在部署包 `dist/`）
- `celatlas_spatial-1.8.0.tar.gz`：源码分发包（包含在部署包 `dist/`）
- `CHECKSUMS.sha256`：部署包内文件校验和

部署包不包含用户 FASTQ、参考索引、样本结果或生产配置。请根据部署指南单独准备这些运行时资产。

## 已知发布边界

- 当前发布包是 Linux x86_64 版本；Windows/macOS 不在本次交付范围。
- v1.8 的完整流程依赖外部 bioinformatics 工具（STAR、featureCounts、samtools、cutadapt）和用户提供的参考索引。
- 合成 FASTQ 仅用于安装、文件发现和 `--dry-run` 冒烟检查，不代表真实生物学结果。
