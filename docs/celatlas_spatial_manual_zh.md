# Celatlas Spatial v1.8.0 中文使用手册

本文档是 Celatlas Spatial v1.8.0 的公开版使用说明，面向 Linux 命令行用户、分析人员和部署人员。内容以当前发布包中的公开 CLI、Python runner 和兼容 shell 入口为准，并参考 v1.7 中文手册整理了输入组织、三种图像模式、两种分析模式和输出结构。

> 文档版本：v1.8.0  
> 适用平台：Linux x86_64  
> 推荐 Python：3.11  
> 说明：命令中的 `/data/...`、样本编号和物种名称都是虚拟示例，请替换成实际路径。

## 1. 软件定位和处理流程

Celatlas Spatial 将空间转录组或单细胞/单核 RNA 数据从 FASTQ 处理到表达矩阵、空间 bin、可选细胞分割、下游聚类和 HTML 报告。它由三层组成：

1. **公开入口层**：`celatlas count`、`celatlas reanalyze`、`celatlas mkreport`。
2. **Python runner 层**：负责配置合并、路径解析、输入预检、步骤计划、恢复和日志。
3. **分析步骤层**：执行 barcode、cutadapt、STAR、featureCounts、UMI count、binSegment、StarDist、分析和报告。

标准步骤如下：

```text
00.sample
  -> 01.barcode
  -> 02.cutadapt
  -> 03.star
  -> 04.featureCounts
  -> 05.count
  -> 06.segment/01.binsegment
  -> 06.segment/02.cellsegment       （可选）
  -> 07.outs/binned_outputs
  -> 07.outs/cellsegmented_outputs   （可选）
  -> 08.report
```

`strna` 保留 barcode 的空间坐标并生成空间 bin；`scrna` 将 barcode 当作独立细胞，跳过空间 bin，直接进行细胞表达分析。`ST`、`SX`、`SN` 表示实验/文库流程；`denovo`、`reanalysis`、`report` 表示从哪里开始执行。

## 2. v1.8 与 v1.7 的主要变化

- 推荐入口统一为 `celatlas`，减少直接拼接 shell 参数造成的路径和恢复问题。
- 增加 `celatlas reanalyze` 和 `celatlas mkreport`，可以在不重复 FASTQ 比对的情况下重做下游分析或只刷新报告。
- 增加 FASTQ manifest。runner 自动记录每个 R1/R2 文件的路径、大小和修改时间，避免把旧 count 结果误用于新 FASTQ。
- 增加 `--dry-run`、`--preflight`、`--show-steps`、`--show-resume-plan` 等检查入口。
- 增加可选 cavity filter、StarDist 细胞分割、DAPI preset 和大图 tiled inference。
- 输出目录统一为 `06.segment`、`07.outs`、`08.report`。旧项目仍可用兼容 shell 入口，但建议新项目使用公开 CLI。
- 发布包不包含内部生产表格、网页调度器、日志、参考索引和实验数据。

## 3. 安装和环境验证

### 3.1 使用正式发布目录

发布包解压目录例如：

```text
/data/software/celatlas-spatial-v1.8.0-linux-x86_64/
```

解压后先校验文件：

```bash
cd /data/software/celatlas-spatial-v1.8.0-linux-x86_64
sha256sum -c CHECKSUMS.sha256
```

使用发布包中的 Conda 安装脚本：

```bash
source /data/miniforge3/etc/profile.d/conda.sh
./release/install.sh
conda activate celatlas18
./release/verify.sh
```

`verify.sh` 应能找到 `celatlas`、`celatlas_spatial`、`STAR`、`featureCounts`、`samtools` 和 `cutadapt`。如果服务器没有 Conda，请先安装 Miniforge/Mambaforge，或按照 `docs/linux_deployment_guide.md` 使用离线环境。

### 3.2 手动激活时的推荐设置

```bash
source /data/miniforge3/etc/profile.d/conda.sh
conda activate celatlas18
export PYTHONNOUSERSITE=1
export MPLCONFIGDIR=/tmp/celatlas_mplconfig
export NUMBA_CACHE_DIR=/tmp/celatlas_numba_cache
mkdir -p "$MPLCONFIGDIR" "$NUMBA_CACHE_DIR"
celatlas --help
celatlas_spatial_runner --help
```

大图或高线程任务还需要足够的临时空间。不要把临时目录放在空间不足的系统分区。

## 4. 发布包目录和外部资产

## 4.1 FASTQ 脱敏 smoke demo

正式发布包包含 `demo_data/fastq/BBV2.4/`，内有 20,000 对配对 reads，以及 `release/create_fastq_demo.py`。该数据用于验证安装、FASTQ 自动发现和 dry-run 计划，不是生物学基准；要完成真实分析仍需完整参考和匹配的空间资产。所有 FASTQ header 均已重写为合成编号；本 demo 的碱基和质量值也已替换为合成内容，仅保留 read 长度和配对结构。manifest 不记录源服务器路径。

在解压后的发布目录执行：

```bash
cd /data/software/celatlas-spatial-v1.8.0-linux-x86_64
source /data/miniforge3/etc/profile.d/conda.sh
conda activate celatlas18
celatlas count \
  --id CELATLAS_DEMO \
  --sample-name Celatlas_Demo \
  --tissue demo \
  --targetdir "$PWD/demo_results/CELATLAS_DEMO" \
  --chemistry BBV2.4 \
  --species Mus_musculus \
  --workflow ST \
  --image gene_expr \
  --mode scrna \
  --fastqs "$PWD/demo_data/fastq/BBV2.4" \
  --reference-dir /data/celatlas/reference \
  --thread 2 \
  --bin 50 \
  --skip-preflight \
  --dry-run
```

真正执行 `count` 时必须提供完整 STAR reference、匹配 chemistry 的输入，以及 `strna` 所需的 barcode/mask。20,000 对 reads 通常不足以产生有意义的比对和 QC。

如需从获得授权的配对 FASTQ 生成新的脱敏片段：

```bash
python release/create_fastq_demo.py \
  --r1 /data/input/sample_R1.fq.gz \
  --r2 /data/input/sample_R2.fq.gz \
  --outdir demo_data/fastq/custom \
  --sample CELATLAS_DEMO \
  --reads 20000
```

工具保留配对 read 的长度并校验 R1/R2 同步，将 header 重写为 `CELATLAS_DEMO:NNNNNNNN/{1,2}`，默认替换碱基和质量值，同时生成不含源路径和原始测序仪信息的 SHA256 清单。只有在获得生物序列再发布授权时才使用 `--sequence-mode copy`。



```text
celatlas-spatial-v1.8.0-linux-x86_64/
├── dist/                 wheel 和 source archive
├── docs/                 中英文手册、部署说明和发布总结
├── envs/                 Conda/Pip 环境定义
├── config-templates/     通用配置模板
├── models/               Swin 模型和可选 StarDist H&E 模型
├── release/              install.sh、verify.sh、运行时资产说明
├── CHECKSUMS.sha256
└── RELEASE_NOTES.md
```

以下内容必须由使用者另外准备：

- FASTQ 文件；
- 空间 barcode 文件或 mask；
- H&E、ssDNA 或 DAPI 图像；
- ST/SN 物种 STAR reference，或 SX panel reference；
- Swin `swin_tiny.pth`（复制到配置的 `src_dir`）；
- 使用细胞分割时的 StarDist 模型。发布包通常带 H&E 模型；`2D_versatile_fluo` 需要按部署说明准备。

参考索引不是只有 FASTA/GTF。目录中还应包含 STAR 索引文件和 Celatlas 使用的配置文件，例如 `Genome`、`SA`、`SAindex` 和 `celatlas_spatial_genome.config`。

### 4.1 从 FASTA/GTF 构建参考

如果手头是 GFF3，先转换为包含 `gene_id` 的 GTF：

```bash
python scripts/normalize_gff3_to_gtf.py \
  /data/reference/raw/genes.gff3 \
  /data/reference/raw/genes.celatlas.gtf
```

再使用 RNA reference builder：

```bash
celatlas_spatial rna mkref \
  --genome_name Mus_musculus_demo \
  --fasta /data/reference/raw/Mus_musculus.fa \
  --gtf /data/reference/raw/genes.celatlas.gtf \
  --outdir /data/celatlas/reference/Mus_musculus \
  --thread 16
```

需要线粒体基因列表时增加 `--mt_gene_list /data/reference/raw/mt_genes.txt`。先检查命令而不构建索引：

```bash
celatlas_spatial rna mkref ... --dry_run
```

`--genome_name`、FASTA、GTF 和输出目录应写入项目记录。SX 的 panel reference 应直接生成到配置中的 `sx_reference_dir`，不要把普通物种 reference 当成 panel reference 使用。

## 5. 三个选择维度

### 5.1 workflow：ST、SX、SN

| workflow | 含义 | 典型用途 | 参考目录 |
|---|---|---|---|
| `ST` | 常规空间转录组 | fresh/frozen 等常规文库 | `reference/<species>` |
| `SX` | FFPE/靶向 panel | FFPE、探针或 targeted panel | `sx_reference_dir` 直接指向 panel reference |
| `SN` | 随机引物方案 | 适用于随机引物建库的空间/单核类流程 | `reference/<species>` |

`FF`、`FFPE` 等名称仍可在兼容 shell 或 Python runner 中作为别名使用，但公开 `celatlas` CLI 只接受 `ST`、`SX`、`SN`。

`SN` 是建库方案标识，表示随机引物方案；它不应仅按“single nucleus（单核）”字面理解。单核 RNA、植物单核或其他采用随机引物建库的样本，均应以实际试剂和建库说明确认后选择 `SN`。

### 5.2 image：gene_expr、ssDNA、HE

| image | 输入 | 组织检测 | 适用场景 |
|---|---|---|---|
| `gene_expr` | 不需要组织图像 | 根据 GEM/UMI 分布计算 | 没有图像时的默认方案 |
| `ssDNA` | `<id>.tif` 等荧光图像 | 传统荧光 mask | ssDNA 或荧光组织成像 |
| `HE` | `<id>_he.tif` 等 H&E 图像 | 基因表达组织检测 + H&E 配准 | 病理和组织形态分析 |

H&E 图像必须带 `_he` 后缀。`<id>.tif` 是 `ssDNA` 的默认命名，不能把普通 H&E 图像误命名成这个形式。

### 5.3 mode：strna、scrna

| mode | 是否使用空间坐标 | 输出重点 |
|---|---|---|
| `strna` | 使用 | bin10/20/50/100、空间聚类、空间图 |
| `scrna` | 不使用 | 细胞/条码表达矩阵、PCA/UMAP/marker |

有 `barcodeToPos.h5` 或 `FilterBarcodes.csv` 并且研究组织空间分布时使用 `strna`。只有单细胞/单核表达分析需求时使用 `scrna`。

## 6. 输入文件准备

### 6.1 推荐目录布局

```text
/data/celatlas/
├── fastq/BBV2.4/
│   ├── SX000293_A1_S1_L001_R1_001.fastq.gz
│   ├── SX000293_A1_S1_L001_R2_001.fastq.gz
│   └── ...
├── mask/
│   ├── SX000293_A1.barcodeToPos.h5
│   ├── SX000293_A1_FilterBarcodes.csv
│   └── SX000293_A1_tissue_bbox.csv
├── images/
│   └── SX000293_A1_he.tif
├── reference/
│   └── Homo_sapiens/
├── src/
│   └── swin_tiny.pth
└── results/
```

`mask-dir` 中最关键的公开输入是 `<id>_FilterBarcodes.csv`。空间 barcode 步骤优先使用 `<id>.barcodeToPos.h5`，若该文件不存在则使用 CSV；`<id>_tissue_bbox.csv` 存在时会用于组织画布范围。

### 6.2 FASTQ 命名和自动发现

`--fastq-name` 默认为 `--id`。runner 按以下顺序查找，并且一旦找到一种布局就不再混用其他布局：

1. 多 lane：`<prefix>_S*_L*_R1_*.fastq.gz` 与对应 R2；也支持 `.fq.gz`。
2. fold：`<prefix>_fold1_1.fq.gz`、`<prefix>_fold1_2.fq.gz`，最多 fold1 到 fold5。
3. 单文件：`<prefix>_1.fq.gz`、`<prefix>_2.fq.gz`。

R1/R2 数量必须匹配，文件必须存在、非空且可读。测序公司提供的前缀不同于芯片编号时：

```bash
--id SX000293_A1 --fastq-name Sequencer_Sample_20261001
```

### 6.3 图像命名

- `HE`：`<id>_he.tif`、`.tiff`、`.png`、`.jpg` 或 `.jpeg`。
- `ssDNA`：`<id>.tif`。运行时 staged 目录中也可以使用相同 stem 的荧光图像。
- `gene_expr`：不需要图像。
- DAPI preset：先准备手工 H&E ROI mask，配准后使用全分辨率 registered TIFF 的蓝通道；不能把 hires 预览图当作分割输入。

### 6.4 mask、ROI 和坐标一致性

空间坐标、组织 bbox、图像和 registered 图像必须位于同一坐标系。StarDist 的 `--stardist-labels` 如果是已有 label TIFF，必须已经在 `FilterBarcodes.csv` 对应的 registered 坐标系中；程序不会根据文件名自动推断变换。

### 6.5 运行时 staging 和手工 ROI mask

公开 `count` 和 `reanalyze` 会把空间输入暂存到：

```text
<targetdir>/06.segment/mask/
```

可识别的手工 mask 文件名包括 `manual_mask.png`、`<id>_manual_mask.png`、`<id>_mask_manual.png`、`mask_manual.png`、`he_manual_mask.png`、`<id>_he_manual_mask.png`、`<id>_he_mask_manual.png` 和 `he_mask_manual.png`。将 mask 放在 `--mask-dir` 后，运行时会自动 staging；也可以直接放入目标目录的 `06.segment/mask`，重分析时已有运行时文件会优先保留。DAPI preset 要求手工 H&E ROI 可用于注册，并要求注册后的全分辨率 TIFF 与 mask 坐标一致。

## 7. 安装后的最小验证

先查看入口：

```bash
celatlas --help
celatlas count --help
celatlas reanalyze --help
celatlas mkreport --help
celatlas_spatial_runner run --help
celatlas_spatial_runner batch --help
```

使用真实数据运行前，先做 dry-run：

```bash
celatlas count \
  --id SX000293_A1 \
  --sample-name Demo_Human_Lung \
  --tissue lung \
  --targetdir /data/celatlas/results/demo/SX000293_A1 \
  --chemistry BBV2.4 \
  --species Homo_sapiens \
  --workflow SX \
  --image HE \
  --mode strna \
  --fastqs /data/celatlas/fastq/BBV2.4 \
  --mask-dir /data/celatlas/mask \
  --image-dir /data/celatlas/images \
  --reference-dir /data/celatlas/reference \
  --config /data/celatlas/config/runner.yaml \
  --dry-run
```

`--targetdir` 必须是绝对路径。dry-run 会打印计划并执行输入检查，不会启动完整分析。

## 8. 公开 CLI

### 8.1 `celatlas count`

从 FASTQ 执行完整流程。必填项：`--id`、`--sample-name`、`--tissue`、`--targetdir`、`--chemistry`、`--species`、`--workflow`、`--image`、`--fastqs`。

```bash
celatlas count \
  --id DEMO_ST000001_A1 \
  --sample-name Demo_Mouse_Brain \
  --tissue brain \
  --targetdir /data/celatlas/results/demo/DEMO_ST000001_A1 \
  --chemistry BBV2.4 \
  --species Mus_musculus \
  --workflow ST \
  --image gene_expr \
  --mode strna \
  --fastqs /data/celatlas/fastq/BBV2.4 \
  --mask-dir /data/celatlas/mask \
  --reference-dir /data/celatlas/reference \
  --thread 32 \
  --bin 50
```

### 8.2 `celatlas reanalyze`

从已有 `05.count` 开始重做 binSegment、可选 post-binSegment、`07.outs` 和报告。它要求 `<targetdir>/05.count/<id>_count_detail.txt` 已存在，不会重新读取 FASTQ 做 STAR 比对。

```bash
celatlas reanalyze \
  --id DEMO_ST000001_A1 \
  --sample-name Demo_Mouse_Brain \
  --tissue brain \
  --targetdir /data/celatlas/results/demo/DEMO_ST000001_A1 \
  --chemistry BBV2.4 \
  --species Mus_musculus \
  --workflow ST \
  --image HE \
  --mode strna \
  --mask-dir /data/celatlas/mask \
  --image-dir /data/celatlas/images \
  --reference-dir /data/celatlas/reference \
  --bin 50
```

`--skip-binsegment` 保留现有 `06.segment/01.binsegment`，从后续分析开始；`--skip-analysis` 保留 `07.outs`，只刷新报告。若要启用 cavity 或细胞分割，reanalysis 会自动从 `06` 规划相应步骤。

### 8.3 `celatlas mkreport`

只生成报告，要求已有完整的下游输出：

```bash
celatlas mkreport \
  --id DEMO_ST000001_A1 \
  --sample-name Demo_Mouse_Brain \
  --tissue brain \
  --targetdir /data/celatlas/results/demo/DEMO_ST000001_A1 \
  --chemistry BBV2.4 \
  --species Mus_musculus \
  --workflow ST \
  --image gene_expr \
  --mode strna
```

## 9. 主要参数说明

### 9.1 路径和元数据

| 参数 | 说明 |
|---|---|
| `--id` | 芯片/样本运行 ID，也是默认 FASTQ 前缀和输出文件 stem |
| `--sample-name` | 生物样本名，写入 metadata 和报告 |
| `--tissue` | 组织名，写入报告 |
| `--targetdir` | 单样本最终目录，必须为绝对路径 |
| `--project` | 项目/案例名；未指定时使用 targetdir 的父目录名 |
| `--fastqs` | FASTQ 目录；`count` 必须提供 |
| `--fastq-name` | FASTQ 实际文件前缀，默认等于 `--id` |
| `--mask-dir` | barcode 和组织坐标目录 |
| `--image-dir` | 原始 H&E/荧光图像目录 |
| `--reference-dir` | ST/SN 参考根目录，SX 时通常指 panel reference |
| `--config` | YAML、JSON 或简单 env 配置 |
| `--profile` | profile 名称或 profile 文件 |
| `--log-dir` | runner 日志和摘要目录 |

### 9.2 线程、分辨率和下游分析

| 参数 | 作用 | 建议 |
|---|---|---|
| `--thread` | 通用 CPU 线程 | 共享服务器按资源保守设置 |
| `--star-thread` | 单独限制 STAR 线程 | 大内存服务器也建议单独设置 |
| `--featurecounts-thread` | featureCounts 线程 | SX 常用 8 |
| `--bin` | 空间 bin，常用 10/20/50/100 | 50 是日常平衡值 |
| `--pixel-size` | 图像像素对应的微米数 | 必须与成像系统匹配 |
| `--cluster-resolution` | 空间 bin Leiden 分辨率 | 越高通常得到更多 cluster |
| `--cell-cluster-resolution` | 细胞 Leiden 分辨率 | 只在细胞分析时使用 |
| `--gem-bin-size` | GEM 热图聚合尺寸 | 10 更细，20 更省资源 |
| `--cell-num` | count 阶段预期细胞/spot 数 | 根据组织规模调整 |
| `--cell-min-genes` / `--cell-min-counts` | 细胞下限过滤 | 低质量数据可适当提高/降低 |
| `--cell-max-genes` / `--cell-max-counts` / `--cell-max-mt` | 上限和线粒体过滤 | 先查看报告再调参 |

bin10 接近单细胞但数据量大，bin20 适合精细结构，bin50 是默认平衡设置，bin100 适合全局模式。`bin` 不等同于细胞直径；它是空间聚合网格的边长。

### 9.3 ssDNA 参数

| 参数 | 作用 |
|---|---|
| `--ssdna-threshold-scale` | 缩放 Otsu 阈值；小于 1 可保留较暗组织边缘 |
| `--ssdna-mask-expand-pixels` | 在 registered 图像上扩张最终组织 mask |
| `--ssdna-min-hole-area` | 只填充小于该面积的 mask 内孔 |
| `--fluorescence-background` | H&E ROI 配准时保留黑色荧光背景 |

这些参数单位是 registered 图像像素或图像阈值比例，不是微米。修改前应保存 dry-run 和报告中的 mask QC。

### 9.4 cavity filter

`--cavity-mode off|qc|apply` 只影响组织 mask 和 bin 输出：

- `off`：不运行。
- `qc`：生成 cavity 过滤结果和 QC，不替换原始 `square_bin`。
- `apply`：在备份原输出后替换用于后续分析的 bin 目录。

相关参数：

```text
--cavity-min-component-ratio   保留的组织连通域相对最大连通域的最小面积比例
--cavity-min-hole-area         填充的小孔面积阈值（全分辨率像素）
--cavity-close-radius          组织 mask 闭运算半径（全分辨率像素）
--preserve-gem-support         在封闭 H&E 空洞内恢复 GEM 支持
--cavity-force                 允许覆盖已有 cavity backup
```

首次建议使用 `qc`，检查 QC 图和统计后再使用 `apply`。`apply` 会保留备份，但仍应先复制样本目录或使用版本化存储。

### 9.5 StarDist 细胞分割

启用方式：

```bash
--enable-cell-segmentation
```

HE + 手工 ROI + DAPI/蓝通道的推荐方式：

```bash
--cell-segmentation-preset dapi
```

`dapi` preset 只能用于 `--image HE --mode strna`。它会要求手工 H&E ROI、启用黑背景、在全分辨率 registered TIFF 上进行 tiled fluorescence segmentation，并运行细胞矩阵和报告步骤。默认通道是 `blue`，也可以显式指定：

```bash
--stardist-fluorescence-channel blue
```

以下高级参数通过兼容 shell、配置/profile 或底层 `celatlas_spatial rna stardistCellSegment` 使用；公共 `celatlas` 命令直接支持的是 `--enable-cell-segmentation`、`--cell-segmentation-preset dapi` 和 `--stardist-fluorescence-channel`。DAPI preset 会自动设置全分辨率 tiled 推理相关默认值。

常用高级参数：

| 参数 | 作用 |
|---|---|
| `--stardist-model` | 模型名，HE 默认 `2D_versatile_he`，荧光/DAPI 默认 `2D_versatile_fluo` |
| `--stardist-model-dir` | 本地模型目录 |
| `--stardist-labels` | 已完成且坐标一致的 label TIFF；提供后跳过推理 |
| `--stardist-image` | 指定分割图像 |
| `--stardist-tissue-bbox` | 指定 label 画布 bbox |
| `--stardist-prob-thresh` | 概率阈值，默认约 0.30 |
| `--stardist-nms-thresh` | NMS 阈值 |
| `--stardist-max-dim` | 普通推理的长边限制；DAPI preset 默认使用原分辨率 |
| `--stardist-scale` | StarDist 内部缩放 |
| `--stardist-n-tiles` | 普通分块网格，如 `4,4` |
| `--stardist-expand-pixels` | 条码归属前扩张 label 的像素数 |
| `--stardist-min-umi` / `--stardist-min-genes` | 保留细胞的最低表达门槛 |
| `--stardist-skip-counts` | 只生成 label/归属，不聚合 count |
| `--stardist-no-label-output` | 不写 slide 尺寸的 label TIFF |
| `--stardist-tiled-inference` | 大图使用可恢复的重叠 tile 推理 |
| `--stardist-tile-size` | tile 边长，默认 4096 |
| `--stardist-tile-overlap` | tile 上下文 halo，默认 256 |
| `--stardist-tile-merge-overlap` | tile label 合并重叠阈值 |
| `--stardist-tile-cache-dir` | tile 缓存目录 |
| `--stardist-exclude-rectangle` | 排除区域 `x0,y0,x1,y1`，可用分号指定多个 |

兼容 shell 的普通 HE 或荧光大图示例：

```bash
bash Celatlas_reanalysis.sh DEMO_HE0002_A1 demo_case BBV2.4 Homo_sapiens HE strna \
  --sample Demo_HE --sampledir /data/celatlas/results/demo/DEMO_HE0002_A1 \
  --image_dir /data/celatlas/images --mask_dir /data/celatlas/mask \
  --reference_dir /data/celatlas/reference \
  --enable-stardist-cell-segment --stardist-tiled-inference \
  --stardist-tile-size 4096 --stardist-tile-overlap 256
```

### 9.6 workflow 的高级参数

这些参数主要通过 `celatlas_spatial_runner run`、CSV/profile 或兼容 shell 入口使用：

| 参数/配置键 | 作用 |
|---|---|
| `--insert-r2` / `--insertR2` | R2 insert 长度；应与测序设计匹配 |
| `--star-match-min` | STAR `outFilterMatchNmin` |
| `--star-match-ratio` | STAR 最小匹配长度比例 |
| `--star-score-ratio` | STAR 最低比对分数比例 |
| `--star-multimap` | STAR `outFilterMultimapNmax` |
| `--bbv4-strna-barcode-mismatch` | BBV4/BBV4_L9 空间 whitelist 允许的 mismatch |
| `--resolve-multigene-umi` | 多基因归属 UMI 按最高 read 支持解析，平票丢弃 |
| `--fastq-root` | 按 chemistry 自动解析 FASTQ 根目录 |
| `--workspace` / `--results-root` | 默认工作区和结果根目录 |
| `--genome-dir` / `--src-dir` | 显式 reference 和模型目录 |
| `--stage-inputs` | 仅执行空间输入 staging 并显示动作 |

SN 的随机引物方案参数通常放在 profile 或 YAML 的 `sn` 部分，例如 `sn_barcode_lownum`、`sn_cutadapt_min_length`、`sn_star_multimap`、`sn_featurecounts_param` 和 `sn_resolve_multigene_umi`。BBV4/HD 灵敏度调试可以使用 `configs/profiles/bbv4_hd_loose.env`，但改变 STAR 过滤阈值前必须保留对照运行和 QC。

## 10. Python runner 高级入口

### 10.1 单样本 step 控制

`celatlas_spatial_runner run` 支持 `shell` 和 `python` 两种 engine。公开 CLI 默认通过 Python runner 规划步骤并调用兼容后端；直接使用 runner 时，`--engine python` 适合查看和执行 Python step plan。

```bash
celatlas_spatial_runner run \
  --engine python \
  --pipeline reanalysis \
  --workflow ST \
  --chip-number DEMO_ST000006_A1 \
  --casno demo_project \
  --sample-name Demo \
  --tissue brain \
  --chemistry BBV2.4 \
  --species Mus_musculus \
  --method gene_expr \
  --mode strna \
  --sampledir /data/celatlas/results/demo/DEMO_ST000006_A1 \
  --from-step 06 \
  --to-step 08.report \
  --show-steps \
  --dry-run
```

常用 step 选择器：

```text
00.sample
01.barcode
02.cutadapt
03.star
04.featureCounts
05.count
06
06.segment/01.binsegment
06.segment/02.cellsegment
07
08.report
```

检查而不执行：

```bash
--show-steps --show-step-status --show-resume-plan --show-execution-preflight --preflight --dry-run
```

`--allow-existing-step-outputs` 只应在确认输出与当前输入相符时使用；否则应使用 reanalysis、manifest 或明确的 reset 选项。

### 10.2 CSV 批处理

批处理接口是公开 Python runner 的通用能力，样本表至少需要：

```text
workflow,chip_number,casno,chemistry,species,method,mode
```

常用字段包括 `run_id`、`enabled`、`sample_name`、`tissue`、`thread`、`bin`、`fastq_dir`、`fastq_name`、`mask_dir`、`image_dir`、`reference_dir`、`resume_existing`、`gene_mask_filter`、`enable_cavity_filter`、`cavity_apply`、`enable_stardist_cell_segment`、`extra_args`。

示例 `samples.csv`：

```csv
run_id,enabled,workflow,sample_name,tissue,chip_number,casno,chemistry,species,method,mode,thread,bin,fastq_dir,mask_dir,image_dir,reference_dir
demo_st,1,ST,Mouse brain,brain,DEMO_ST0007_A1,demo_case,BBV2.4,Mus_musculus,gene_expr,strna,32,50,/data/celatlas/fastq/BBV2.4,/data/celatlas/mask,/data/celatlas/images,/data/celatlas/reference
demo_sx,1,SX,Human lung,lung,DEMO_SX0008_A1,demo_case,BBV2.4,Homo_sapiens,HE,strna,32,50,/data/celatlas/fastq/BBV2.4,/data/celatlas/mask,/data/celatlas/images,/data/celatlas/reference
```

运行：

```bash
celatlas_spatial_runner batch \
  --samples /data/celatlas/config/samples.csv \
  --config /data/celatlas/config/runner.yaml \
  --preflight --dry-run --show-steps
```

只运行指定行：

```bash
--only demo_sx
```

批处理失败时可用 `--stop-on-error` 在首个失败样本停止。

### 10.3 低层 cavity 和 StarDist 入口

需要单独重做 post-binSegment 时，可以直接调用 RNA 子命令：

```bash
celatlas_spatial rna cavityFilter \
  --sampledir /data/celatlas/results/demo/DEMO_HE0002_A1 \
  --bins 20,50 --mask-preset default --dry-run
```

已有同一 registered 坐标系的 label 时，可以跳过推理：

```bash
celatlas_spatial rna stardistCellSegment \
  --sampledir /data/celatlas/results/demo/DEMO_HE0002_A1 \
  --labels /data/celatlas/labels/DEMO_HE0002_A1.labels.tif \
  --skip-inference --min-umi 10 --min-genes 5
```

大图从原分辨率推理：

```bash
celatlas_spatial rna stardistCellSegment \
  --sampledir /data/celatlas/results/demo/DEMO_HE0002_A1 \
  --tiled-inference --tile-size 4096 --tile-overlap 256 \
  --fluorescence-channel blue --tile-cache-dir /data/celatlas/tile_cache
```

低层入口仍要求 `FilterBarcodes.csv`、`tissue_bbox.csv`、registered TIFF 和 label canvas 坐标一致。

## 11. 配置文件、profile 和优先级

### 11.1 YAML 配置

```yaml
paths:
  workspace: /data/celatlas
  results_root: /data/celatlas/results
  fastq_root: /data/celatlas/fastq
  mask_dir: /data/celatlas/mask
  image_dir: /data/celatlas/images
  reference_dir: /data/celatlas/reference
  sx_reference_dir: /data/celatlas/reference/Homo_sapiens_wtpanel
  src_dir: /data/celatlas/src

defaults:
  default_thread: 32
  default_star_thread: 12
  default_bin: 50
  default_pixel_size: 0.5
  gem_bin_size: 20

sx:
  featurecounts_thread: 8

env:
  CELATLAS_ENV_NAME: celatlas18
  CELATLAS_AUTO_ACTIVATE: "1"
  CELATLAS_USE_ENV_PATH: "1"
```

### 11.2 env 配置

```bash
export CELATLAS_WORKSPACE=/data/celatlas
export CELATLAS_REFERENCE_DIR=/data/celatlas/reference
export CELATLAS_MASK_DIR=/data/celatlas/mask
export CELATLAS_IMAGE_DIR=/data/celatlas/images
export CELATLAS_FASTQ_ROOT=/data/celatlas/fastq
export CELATLAS_RESULTS_ROOT=/data/celatlas/results
export CELATLAS_SRC_DIR=/data/celatlas/src
```

用 `--config /path/to/site.env` 传入。解析器读取 `KEY=value` 或 `export KEY=value`，不会执行配置文件中的 shell 控制流。

### 11.3 合并规则

通常优先级从低到高为：站点配置 → profile → workflow/method 配置 → CSV 行或命令行值。命令行显式值应作为最终确认值。配置中的 `paths`、`defaults`、`env`、`methods`、`ST/SX/SN` 和 profile 可分别保存通用路径、默认值和流程专用设置。

## 12. 恢复、重跑和 FASTQ manifest

### 12.1 安全恢复

首次运行建议：

```bash
celatlas count ... --dry-run
celatlas count ... --resume-existing
```

manifest 成功提交后位于：

```text
<targetdir>/.fastq_inputs.tsv
```

任务执行期间尝试的 manifest 位于 `.fastq_inputs.current.tsv`。只有成功运行后 current 才会提升为 baseline，因此失败重跑不会把失败输入伪装成已验证输入。

### 12.2 FASTQ 改变后的重跑

`celatlas count` 会自动建立和检查 FASTQ manifest。检测到 FASTQ 增加、顺序变化或文件大小/时间变化时，runner 会阻止直接复用依赖 FASTQ 的旧输出。需要明确重置：

```bash
celatlas count ... --reset-fastq-outputs
```

希望保留旧结果时：

```bash
celatlas count ... --reset-fastq-outputs --archive-fastq-outputs
```

旧目录会被放入样本目录下的 `.rerun_archive`。不要在未确认 manifest 的情况下使用 `--allow-existing-step-outputs`。

### 12.3 下游重分析

图像、mask、bin、cavity 或细胞分割发生变化时优先使用 `reanalyze`，避免重复 STAR。对于只修改报告模板或报告输入的情况使用 `mkreport`。

## 13. 输出目录和结果解释

```text
<targetdir>/
├── 00.sample/                         样本信息和初始统计
├── 01.barcode/                        barcode 提取结果
├── 02.cutadapt/                       接头去除和长度过滤
├── 03.star/                           STAR BAM 和比对日志
├── 04.featureCounts/                  基因计数 BAM/表格
├── 05.count/                          UMI count、raw/filtered 矩阵
│   ├── <id>_count_detail.txt
│   └── <id>_filtered_feature_bc_matrix/
├── 06.segment/
│   ├── 01.binsegment/
│   │   ├── images/                    mask、registered image、预览图
│   │   └── square_bin/                bin10/20/50/100/Raw
│   └── 02.cellsegment/                可选 StarDist label 和归属
├── 07.outs/
│   ├── binned_outputs/                空间降维、cluster、图像
│   └── cellsegmented_outputs/         可选细胞级分析
├── 08.report/                         报告步骤日志/中间文件
├── <id>_spatial_analysis_report.html
├── .fastq_inputs.tsv                  成功运行的 FASTQ 基线
└── pipeline.log / runner 日志
```

空间矩阵通常采用 10X 风格稀疏格式：

```text
filtered_feature_bc_matrix/
├── barcodes.tsv.gz
├── features.tsv.gz
└── matrix.mtx.gz
```

可用 Seurat 的 `Read10X()` 或 Scanpy 的 `sc.read_10x_mtx()` 读取。报告中的 UMI、基因数、比对率、细胞数、空间聚类和图像叠加应结合实验设计解释；低质量 barcode、组织外区域和图像配准误差会直接影响下游结果。

## 14. 八个虚拟运行 demo

### Demo 1：gene_expr 完整流程

```bash
celatlas count \
  --id DEMO_ST000001_A1 --sample-name Demo_Mouse_Brain --tissue brain \
  --targetdir /data/celatlas/results/demo/DEMO_ST000001_A1 \
  --chemistry BBV2.4 --species Mus_musculus --workflow ST \
  --image gene_expr --mode strna \
  --fastqs /data/celatlas/fastq/BBV2.4 \
  --mask-dir /data/celatlas/mask --reference-dir /data/celatlas/reference \
  --thread 32 --bin 50 --dry-run
```

### Demo 2：SX + H&E + StarDist

```bash
celatlas count \
  --id DEMO_SX000002_A1 --sample-name Demo_Human_Lung --tissue lung \
  --targetdir /data/celatlas/results/demo/DEMO_SX000002_A1 \
  --chemistry BBV2.4 --species Homo_sapiens --workflow SX \
  --image HE --mode strna \
  --fastqs /data/celatlas/fastq/BBV2.4 \
  --mask-dir /data/celatlas/mask --image-dir /data/celatlas/images \
  --reference-dir /data/celatlas/reference \
  --enable-cell-segmentation --dry-run
```

### Demo 3：DAPI preset

```bash
celatlas reanalyze \
  --id DEMO_ST000003_A1 --sample-name Demo_DAPI_Liver --tissue liver \
  --targetdir /data/celatlas/results/demo/DEMO_ST000003_A1 \
  --chemistry BBV2.4 --species Mus_musculus --workflow ST \
  --image HE --mode strna --mask-dir /data/celatlas/mask \
  --image-dir /data/celatlas/images --reference-dir /data/celatlas/reference \
  --cell-segmentation-preset dapi --stardist-fluorescence-channel blue \
  --dry-run
```

### Demo 4：ssDNA cavity QC

```bash
celatlas reanalyze \
  --id DEMO_ST000004_A1 --sample-name Demo_Mouse_Kidney --tissue kidney \
  --targetdir /data/celatlas/results/demo/DEMO_ST000004_A1 \
  --chemistry BBV2.4 --species Mus_musculus --workflow ST \
  --image ssDNA --mode strna --mask-dir /data/celatlas/mask \
  --image-dir /data/celatlas/images --reference-dir /data/celatlas/reference \
  --cavity-mode qc --ssdna-threshold-scale 0.85 \
  --ssdna-mask-expand-pixels 16 --dry-run
```

### Demo 5：只生成报告

```bash
celatlas mkreport \
  --id DEMO_SX000005_A1 --sample-name Demo_Report --tissue colon \
  --targetdir /data/celatlas/results/demo/DEMO_SX000005_A1 \
  --chemistry BBV2.4 --species Homo_sapiens --workflow SX \
  --image HE --mode strna
```

### Demo 6：Python runner 的 step 控制

```bash
celatlas_spatial_runner run --engine python --pipeline reanalysis \
  --workflow ST --chip-number DEMO_ST000006_A1 --casno demo_project \
  --sample-name Demo --tissue brain --chemistry BBV2.4 \
  --species Mus_musculus --method gene_expr --mode strna \
  --sampledir /data/celatlas/results/demo/DEMO_ST000006_A1 \
  --from-step 06 --to-step 08.report --show-steps --dry-run
```

### Demo 7：配置文件

```bash
celatlas count \
  --id DEMO_CONFIG_001 --sample-name Demo_Config --tissue brain \
  --targetdir /data/celatlas/results/demo/DEMO_CONFIG_001 \
  --chemistry BBV2.4 --species Mus_musculus --workflow ST \
  --image gene_expr --fastqs /data/celatlas/fastq/BBV2.4 \
  --config /data/celatlas/config/runner.yaml --dry-run
```

### Demo 8：恢复和 FASTQ manifest

```bash
celatlas count \
  --id DEMO_RESUME_001 --sample-name Demo_Resume --tissue brain \
  --targetdir /data/celatlas/results/demo/DEMO_RESUME_001 \
  --chemistry BBV2.4 --species Mus_musculus --workflow ST \
  --image gene_expr --fastqs /data/celatlas/fastq/BBV2.4 \
  --resume-existing --dry-run
```

## 15. 兼容 shell 入口

旧脚本仍随源码和发布包提供：

```text
Celatlas.sh              ST 常规流程
Celatlas_FFPE.sh         SX/FFPE 流程
Celatlas_SN.sh           SN 流程
Celatlas_reanalysis.sh   下游重分析
Celatlas_run.sh          配置/调度兼容入口
```

旧式位置参数示例：

```bash
bash Celatlas.sh DEMO ST110001_A1 project01 BBV2.4 Mus_musculus HE strna \
  --thread 32 --bin 50
```

新项目建议使用 `celatlas count`，因为公开 CLI 会校验必填 metadata、绝对输出路径、FASTQ 和输入预检。shell 入口适合已有调度系统迁移；脚本的 `--reference_dir`、`--mask_dir`、`--image_dir`、`--fastq_dir` 和 `--fastq_name` 仍可用于共享存储布局。

## 16. 常见问题和排错

### Q1：提示找不到 FASTQ

确认 `--fastqs` 是目录而不是单个文件，并确认 `--fastq-name` 与前缀一致：

```bash
find /data/celatlas/fastq/BBV2.4 -maxdepth 1 -type f | sort | head
```

注意 R1/R2 必须成对；多 lane 文件不要混入其他样本前缀。

### Q2：提示找不到 H&E

确认文件名是 `<id>_he.tif`，并且 `--image-dir` 指向该目录。`.tiff`、`.png`、`.jpg`、`.jpeg` 也支持。

### Q3：为什么 `strna` 需要 mask，而 `scrna` 不需要？

`strna` 要把 barcode 映射到空间画布，因此需要 `FilterBarcodes.csv` 或 `barcodeToPos.h5`。`scrna` 直接分析表达矩阵，不执行空间 bin。

### Q4：`reanalyze` 提示缺少 count_detail

确认 `<targetdir>/05.count/<id>_count_detail.txt` 存在。若只有 FASTQ，应使用 `count`；若 count 输出在其他位置，修改 `--targetdir` 或配置中的 sampledir。

### Q5：任务中断后应该怎么继续？

先运行带 `--dry-run --show-resume-plan --preflight` 的同一命令，检查输入和计划。FASTQ 未变化时可使用 `--resume-existing`；若输入发生变化，先用 `--reset-fastq-outputs`，需要留档时再加 `--archive-fastq-outputs`。

### Q6：DAPI 分割结果为空或错位

确认使用 `--cell-segmentation-preset dapi`，已有 registered TIFF 和手工 H&E ROI，`FilterBarcodes.csv` 与 label 在同一 registered 坐标系，并检查蓝通道是否正确。不要把低分辨率 `tissue_hires_image.png` 当作 DAPI 推理输入。

### Q7：StarDist 大图内存不足

减少 `--thread`，使用 `--stardist-tiled-inference`，设置合理的 `--stardist-tile-size` 和缓存目录，并确保缓存所在磁盘有空间。可先用 `--stardist-no-label-output` 降低输出占用，但这不会减少推理本身的显存/内存需求。

### Q8：cavity apply 是否会覆盖原始结果？

`apply` 会在替换前创建 backup；仍建议先用 `--cavity-mode qc` 检查 QC，再复制样本目录或使用归档存储。已有 backup 时如确实要重新生成才使用 `--cavity-force`。

### Q9：报告中的空间图和组织图不重合

检查 `pixel-size`、barcode 坐标、tissue bbox、图像方向和 registered metadata。不要仅根据文件名判断坐标是否一致；重新配准后应通过 `reanalyze` 重建后续输出。

### Q10：如何收集故障信息？

保留以下信息：

```bash
celatlas --help
python --version
celatlas count ... --dry-run --show-steps --preflight
```

同时提供样本目录中的 `pipeline.log`、runner 日志、失败报告、软件版本、参考索引版本和输入文件 manifest。不要上传 FASTQ、患者信息或内部路径到公开 issue。

## 17. 性能、资源和复现建议

- STAR 的内存主要由 reference 和排序阶段决定；`--star-thread` 不宜盲目等于机器总核数。
- bin10、tiled inference 和全分辨率 label 会明显增加内存、临时空间和 I/O。
- 共享服务器应限制 `--thread`，并为 `MPLCONFIGDIR`、`NUMBA_CACHE_DIR` 和 tile cache 指定可写目录。
- 每个项目保存 release archive、`CHECKSUMS.sha256`、配置文件、参考索引版本、模型版本、命令行和 FASTQ manifest。
- 将参考、图像、FASTQ 和结果放在发布目录之外，便于升级和回滚。
- 修改阈值后应保留 dry-run、QC 图、报告和参数记录，避免只保留最终矩阵。

## 18. 已知限制

- 公开发行包不提供物种或 panel 的 STAR reference；用户必须准备与文库 chemistry 匹配的索引。
- 细胞分割质量取决于图像通道、配准、模型和组织类型，StarDist 输出需要人工 QC。
- `gene_expr` 的组织范围由表达信号推断，不能替代真实组织图像验证。
- `strna` 的空间图和 `scrna` 的 UMAP/cluster 是分析结果，不代表自动完成细胞类型注释；marker 和生物学解释需由用户完成。
- 旧 v1.7 文档中的内部生产命令、旧目录名称和未公开参数不属于 v1.8 开源接口。

## 19. 相关文件

- [English user manual](celatlas_spatial_manual_en.md)
- [Linux deployment guide](linux_deployment_guide.md)
- [v1.8 release summary](celatlas-v1.8-release-summary.md)
- [runner config template](../configs/runner.yaml.example)
- [environment definition](../envs/celatlas18.yml)
