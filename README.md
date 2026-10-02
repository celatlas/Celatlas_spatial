# Celatlas Spatial v1.8.0

[![Version](https://img.shields.io/badge/version-1.8.0-1f883d.svg)](https://github.com/celatlas/Celatlas_spatial/releases/tag/v1.8.0)
[![License](https://img.shields.io/badge/license-MIT-2ea44f.svg)](LICENSE.txt)
[![Platform](https://img.shields.io/badge/platform-Linux%20x86__64-555.svg)](#安装与部署)

Celatlas Spatial 是赛图生信面向空间转录组与相关 RNA 流程的 Linux 分析管线。它把 FASTQ、空间坐标和可选的 HE/ssDNA 图像组织成可复现的计数矩阵、空间 bin、细胞级结果和 HTML 报告。

**v1.8.0 是独立的重大更新版本。** 它保留 v1.7 的兼容 shell 入口，同时引入统一的 `celatlas` CLI、ST/SX/SN 工作流、预检与断点续跑、重新分析/报告生成，以及基于 StarDist 的细胞分割和 cell-level 分析。v1.8 的源码位于 `v1.8.0` 分支，发布包位于 [GitHub Releases](https://github.com/celatlas/Celatlas_spatial/releases/tag/v1.8.0)。

## v1.8.0 主要更新

- **统一公共 CLI**：`celatlas count`、`celatlas reanalyze`、`celatlas mkreport`，并保留 `celatlas_spatial` 高级入口。
- **三类工作流**：支持 ST、SX、SN 文库设计；`strna` 与 `scrna` 分析模式可按实验方案选择。
- **完整图像流程**：gene expression、ssDNA 和 HE 图像输入；组织分割、HE 配准、空间 bin 和交互式报告。
- **细胞级分析**：可选 StarDist H&E/荧光细胞分割，生成 cell matrix、AnnData、空间坐标、QC 和细胞级聚类/marker 结果。
- **可复现运行**：FASTQ 布局预检、运行计划、步骤级执行、失败报告、断点续跑和已有 count 结果的 reanalysis。
- **可交付发布包**：提供 Conda 环境定义、Linux 安装/验证脚本、双语手册和脱敏合成 FASTQ 冒烟数据。

## 快速开始

### 1. 下载发布包

从 [v1.8.0 Release](https://github.com/celatlas/Celatlas_spatial/releases/tag/v1.8.0) 下载：

```text
celatlas-spatial-v1.8.0-linux-x86_64.tar.gz
```

```bash
tar -xzf celatlas-spatial-v1.8.0-linux-x86_64.tar.gz
cd celatlas-spatial-v1.8.0-linux-x86_64
./release/install.sh
conda activate celatlas18
./release/verify.sh
```

安装脚本只安装软件和 Python/Conda 依赖。FASTQ、barcode/mask、参考基因组索引和项目数据需要按部署指南单独准备。

### 2. 先做计划和输入检查

```bash
celatlas count \
  --id SX000293_A1 \
  --sample-name Human_FFPE_001 \
  --tissue lung \
  --targetdir /data/celatlas/results/SXV1.1test/SX000293_A1 \
  --chemistry BBV2.4 \
  --species Homo_sapiens \
  --workflow SX \
  --image HE \
  --mode strna \
  --fastqs /data/celatlas/fastq/BBV2.4 \
  --reference-dir /data/celatlas/reference/Homo_sapiens_wtpanel \
  --dry-run
```

确认计划无误后，移除 `--dry-run` 执行完整流程。需要细胞分割时增加 `--enable-cell-segmentation`；需要 GPU/模型和图像输入时，参见部署指南与中文手册。

### 3. 使用已有结果

```bash
# 从已有 count 结果重新进行下游分析
celatlas reanalyze --targetdir /data/celatlas/results/SXV1.1test/SX000293_A1

# 只重新生成 HTML 报告
celatlas mkreport --targetdir /data/celatlas/results/SXV1.1test/SX000293_A1
```

## 入口与适用场景

| 入口 | 用途 |
|---|---|
| `celatlas count` | 新样本：FASTQ → count → 空间分析/报告 |
| `celatlas reanalyze` | 使用已有 count 结果重新做下游分析 |
| `celatlas mkreport` | 只生成或刷新 HTML 报告 |
| `celatlas_spatial run` | 配置文件/CSV 批处理和步骤级控制 |
| `Celatlas.sh`、`Celatlas_FFPE.sh`、`Celatlas_SN.sh`、`Celatlas_reanalysis.sh` | v1.7 兼容 shell 入口 |

## 输出内容

常规空间流程会生成 QC、组织/图像配准结果、不同分辨率的 square-bin 矩阵、UMAP/聚类/marker 和交互式 HTML 报告。启用细胞分割后，还会生成：

- `cell_matrix/`：10X 兼容的细胞级矩阵
- `*_cell_adata.h5ad`：AnnData 数据
- `*_cell_positions.tsv`、`*_cell_metadata.tsv`：细胞坐标和 QC 元数据
- 细胞边界/overlay、cell-level QC、聚类、空间图和 marker 结果

项目数据、参考索引和运行结果不应提交到源码仓库。

## 安装与部署

- [中文用户手册](docs/celatlas_spatial_manual_zh.md)
- [English user manual](docs/celatlas_spatial_manual_en.md)
- [Linux 部署指南](docs/linux_deployment_guide.md)
- [Conda 环境定义](envs/celatlas18.yml)
- [Pip 依赖清单](envs/celatlas18.requirements.txt)
- [发布包说明](RELEASE_NOTES.md)

最低要求：Linux x86_64、Python 3.11、可用的 Conda/Miniforge 环境。完整流程还需要 STAR、featureCounts、samtools 和 cutadapt；HE/细胞分割流程需要相应图像、模型与计算资源。

## 发布文件

源码仓库用于版本化代码、文档和小型示例；Linux 部署包、模型和校验文件作为 Release 附件发布：

```text
celatlas-spatial-v1.8.0-linux-x86_64.tar.gz
```

下载后请先阅读包内 `RELEASE_NOTES.md`、`release/MISSING_RUNTIME_ASSETS.md`，并保存 `CHECKSUMS.sha256`。

## 许可证与反馈

本项目使用 [MIT License](LICENSE.txt)。问题、功能建议和安装反馈请提交到 [GitHub Issues](https://github.com/celatlas/Celatlas_spatial/issues)。
