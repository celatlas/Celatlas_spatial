#!/bin/bash
set -e

# Input validation and parameter processing
if [ $# -eq 7 ]; then
    # Full usage: sample_name chip_number casno chemistry Species method mode
    sample_name="$1"
    chip_number="$2"
    casno="$3"
    chemistry="$4"
    Species="$5"
    method="$6"
    mode="$7"
    # Use chip_number as the internal sample parameter
    sample="$chip_number"
elif [ $# -eq 6 ]; then
    # Usage without sample_name: chip_number casno chemistry Species method mode
    sample_name=""
    chip_number="$1"
    casno="$2"
    chemistry="$3"
    Species="$4"
    method="$5"
    mode="$6"
    # Use chip_number as the internal sample parameter
    sample="$chip_number"
else
    echo "Usage:"
    echo "  $0 <sample_name> <chip_number> <casno> <chemistry> <Species> <method> <mode>"
    echo "  $0 <chip_number> <casno> <chemistry> <Species> <method> <mode>"
    echo ""
    echo "Parameters:"
    echo "  sample_name  - Sample name (optional in 6-parameter mode)"
    echo "  chip_number  - Chip number"
    echo "  casno        - Case number for result storage (e.g., CAS250801)"
    echo "  chemistry    - Chemistry type"
    echo "  Species      - Species name"
    echo "  method       - Processing method (image/gene_expr)"
    echo "  mode         - Analysis mode"
    exit 1
fi

# Export parameters for child processes
export CELATLAS_SAMPLE_NAME="$sample_name"
export CELATLAS_CHIP_NUMBER="$chip_number"

# Validate method (skip validation for scRNA mode)
if [[ "$mode" != "scrna" && "$method" != "image" && "$method" != "gene_expr" ]]; then
    echo "Error: Invalid method specified. Must be 'image' or 'gene_expr'."
    exit 1
fi

# Activate environment
#conda activate celatlas_spatial

# Define constants
bin=50
pixelSize=0.5
insertR2=150
cell_num=50000
feature_type="gene"
thread=96
tmp="yes"
chemistryPattern="C4L15C4L15C4U10T18"

# Multi-file processing configuration
# Set MAX_PARALLEL_FILES environment variable to override default
# Example: export MAX_PARALLEL_FILES=3 before running this script
MAX_PARALLEL_FILES=${MAX_PARALLEL_FILES:-3}  # Default: process 2 files in parallel

workspace_dir=${CELATLAS_WORKSPACE:-$(readlink -f "/your/path")}
segImageDir="${workspace_dir}/binSegment/${sample}"
image_dir="${workspace_dir}/images"
fastq_dir="${workspace_dir}/fastq/$chemistry"
mask_dir="${workspace_dir}/ST_mask"
reference_dir="${workspace_dir}/reference"
src_dir="${workspace_dir}/src"
rawdata_dir="${workspace_dir}/rawdata/$sample"
sampledir="${workspace_dir}/results/${casno}/${sample}"
genomeDir="${reference_dir}/${Species}"

# Create output directories first
mkdir -p "${rawdata_dir}" "${sampledir}/01.barcode" "${sampledir}"

# Logging
log_file="${sampledir}/pipeline.log"
exec > >(tee -a "$log_file") 2>&1
echo "Starting pipeline at $(date)"

if [[ "$mode" == "scrna" ]]; then
    # scRNA mode does not require spatial position files
    required_files=()
elif [[ "$method" == "image" ]]; then
    required_files=(
        "${image_dir}/${sample}.tif"
        "${mask_dir}/${sample}.barcodeToPos.h5"
        "${mask_dir}/${sample}_FilterBarcodes.csv"
        "${mask_dir}/${sample}_tissue_bbox.csv"
    )
elif [[ "$method" == "gene_expr" ]]; then
    required_files=(
        "${mask_dir}/${sample}.barcodeToPos.h5"
        "${mask_dir}/${sample}_FilterBarcodes.csv"
        "${mask_dir}/${sample}_tissue_bbox.csv"
    )
else
    echo "Error: Unknown method: $method. Valid methods are 'image' and 'gene_expr'."
    exit 1
fi

# Copy required files only if not in scRNA mode
if [[ "$mode" != "scrna" ]]; then
    for file in "${required_files[@]}"; do
        if [[ ! -f "$file" ]]; then
            echo "Error: Required file not found: $file"
            exit 1
        fi
        cp "$file" "${rawdata_dir}"
    done
fi

# Set ulimit
ulimit -n 10240

# Define function: run command and check success
run_command() {
    echo "Running: $@"
    "$@"
    if [ $? -ne 0 ]; then
        echo "Error: Command failed: $@"
        exit 1
    fi
}

# Define function: run barcode step
run_barcode_step() {
    # Check if FASTQ files exist in rawdata_dir
    if [[ -f "${rawdata_dir}/${sample}_1.fq.gz" && -f "${rawdata_dir}/${sample}_2.fq.gz" ]]; then
        echo "FASTQ files already exist in rawdata directory. Copying to barcode directory..."
        cp "${rawdata_dir}/${sample}_1.fq.gz" "${sampledir}/01.barcode/"
        cp "${rawdata_dir}/${sample}_2.fq.gz" "${sampledir}/01.barcode/"
    # Check if FASTQ files exist in barcode directory
    elif [[ -f "${sampledir}/01.barcode/${sample}_1.fq.gz" && -f "${sampledir}/01.barcode/${sample}_2.fq.gz" ]]; then
        echo "FASTQ files already exist in barcode directory. Skipping barcode step."
    # If neither exists, run the barcode command from fastq_dir
    else
        echo "FASTQ files not found in rawdata or barcode directories. Running barcode step from fastq_dir..."
        
        # Check for fold files (multi-file processing)
        fq1_files=""
        fq2_files=""
        max_parallel=${MAX_PARALLEL_FILES}  # Use configured value
        
        # Look for fold files first
        fold_files_exist=false
        for fold in fold1 fold2 fold3 fold4 fold5; do
            if [[ -f "${fastq_dir}/${sample}_${fold}_1.fq.gz" && -f "${fastq_dir}/${sample}_${fold}_2.fq.gz" ]]; then
                if [[ -n "$fq1_files" ]]; then
                    fq1_files="${fq1_files},"
                    fq2_files="${fq2_files},"
                fi
                fq1_files="${fq1_files}${fastq_dir}/${sample}_${fold}_1.fq.gz"
                fq2_files="${fq2_files}${fastq_dir}/${sample}_${fold}_2.fq.gz"
                fold_files_exist=true
                echo "Found fold file: ${sample}_${fold}_1.fq.gz"
            fi
        done
        
        # If no fold files found, use single files
        if [[ "$fold_files_exist" == false ]]; then
            if [[ -f "${fastq_dir}/${sample}_1.fq.gz" && -f "${fastq_dir}/${sample}_2.fq.gz" ]]; then
                fq1_files="${fastq_dir}/${sample}_1.fq.gz"
                fq2_files="${fastq_dir}/${sample}_2.fq.gz"
                echo "Using single FASTQ files: ${sample}_1.fq.gz, ${sample}_2.fq.gz"
            else
                echo "Error: No FASTQ files found for sample ${sample}"
                exit 1
            fi
        else
            echo "Found multiple fold files, using multi-file parallel processing with max_parallel=${max_parallel}"
        fi
        
        # Set whitelist parameter based on mode
        if [[ "$mode" == "scrna" ]]; then
            whitelist_param=""  # scRNA mode uses chemistry-based whitelist automatically
            echo "scRNA mode: Using chemistry-based whitelist for barcode extraction"
        else
            whitelist_param="--whitelist ${rawdata_dir}/${sample}.barcodeToPos.h5"
        fi

        # Run barcode command with appropriate parameters
        if [[ "$fold_files_exist" == true ]]; then
            # Multi-file processing with parallelization
            run_command celatlas_spatial rna barcode --outdir "${sampledir}/01.barcode" --sample "${sample}" \
                --thread "${thread}" --chemistry "${chemistry}" --pattern "${chemistryPattern}" \
                ${whitelist_param} --mode "${mode}" --lowNum 2 --gzip --output_R1 --resume \
                --max_parallel_files "${max_parallel}" \
                --fq1 "${fq1_files}" --fq2 "${fq2_files}"
        else
            # Single file processing (original behavior)
            run_command celatlas_spatial rna barcode --outdir "${sampledir}/01.barcode" --sample "${sample}" \
                --thread "${thread}" --chemistry "${chemistry}" --pattern "${chemistryPattern}" \
                ${whitelist_param} --mode "${mode}" --lowNum 2 --gzip --output_R1 --resume \
                --fq1 "${fq1_files}" --fq2 "${fq2_files}"
        fi
    fi
}

run_image_branch() {
    echo "Running image branch..."
    run_command celatlas_spatial rna binSegment --outdir ${sampledir}/06.binSegment --sample ${sample} \
        --thread ${thread} --genomeDir ${genomeDir} --pixel-size ${pixelSize} --input ${rawdata_dir} \
        --segment  --tif ${rawdata_dir}/${sample}.tif --bs_out ${segImageDir} --model ${src_dir}/swin_tiny.pth --method ${method} \
        --count --count_detail ${sampledir}/05.count/${sample}_count_detail.txt
}

run_gene_expr_branch() {
    echo "Running gene_expr branch..."
    run_command celatlas_spatial rna binSegment --outdir "${sampledir}/06.binSegment" --sample "${sample}" \
        --thread "${thread}" --genomeDir "${reference_dir}/${Species}" --pixel-size "${pixelSize}" --input "${rawdata_dir}" \
        --segment --model "${src_dir}/swin_tiny.pth" --method "gene_expr" --count \
        --count_detail "${sampledir}/05.count/${sample}_count_detail.txt"
}

# Define function: generate spatial analysis report
generate_spatial_report() {
    echo "Generating spatial analysis report..."
    
    # Check if the report generator exists in the installed package
    if command -v celatlas_spatial_report >/dev/null 2>&1; then
        # Use the installed command
        run_command celatlas_spatial_report \
            "${sampledir}" \
            "${sample}" \
            --chemistry "${chemistry}" \
            --species "${Species}" \
            --output-filename "${sample}_spatial_analysis_report.html" \
            --verbose
        
        echo "Spatial analysis report generated successfully: ${sampledir}/${sample}_spatial_analysis_report.html"
    else
        # Fallback to direct Python execution
        report_generator="/mnt/strna/celatlas_spatial/test/celatlas_spatial/celatlas_spatial/tools/spatial_report_generator.py"
        
        if [[ -f "$report_generator" ]]; then
            # Run the report generator directly
            run_command python3 "$report_generator" \
                "${sampledir}" \
                "${sample}" \
                --chemistry "${chemistry}" \
                --species "${Species}" \
                --output-filename "${sample}_spatial_analysis_report.html" \
                --verbose
            
            echo "Spatial analysis report generated successfully: ${sampledir}/${sample}_spatial_analysis_report.html"
        else
            echo "Warning: Report generator not found. Please install the package or check the path."
            echo "Skipping report generation."
        fi
    fi
}

# Run pipeline steps

# 01.sample
run_command celatlas_spatial rna sample --outdir "${sampledir}/00.sample" --sample "${sample}" \
    --thread "${thread}" --chemistry "${chemistry}" --fq1 "${fastq_dir}/${sample}_1.fq.gz"

# 02.barcode
run_barcode_step

# 03.cutadapt
run_command celatlas_spatial rna cutadapt --outdir "${sampledir}/02.cutadapt" --sample "${sample}" \
    --thread "${thread}" --minimum_length 20 --nextseq_trim 20 --gzip --overlap 10 --insert "${insertR2}" \
    --fq "${sampledir}/01.barcode/${sample}_2.fq.gz"

# 04.STAR
run_command celatlas_spatial rna star --outdir "${sampledir}/03.star" --sample "${sample}" \
    --thread "${thread}" --genomeDir "${reference_dir}/${Species}" \
    --STAR_param '--outFilterMatchNminOverLread 0.2 --outFilterScoreMinOverLread 0.2' \
    --outFilterMultimapNmax 1 --starMem 30 \
    --fq "${sampledir}/02.cutadapt/${sample}_clean_2.fq.gz"

# 05.feature_counts
run_command celatlas_spatial rna featureCounts --outdir "${sampledir}/04.featureCounts" --sample "${sample}" \
    --thread "${thread}" --gtf_type "${feature_type}" --genomeDir "${reference_dir}/${Species}" --featureCounts_param '-s 1' \
    --input "${sampledir}/03.star/${sample}_Aligned.sortedByCoord.out.bam"

# 06.counts
run_command celatlas_spatial rna count --outdir "${sampledir}/05.count" --sample "${sample}" \
    --thread "${thread}" --genomeDir "${reference_dir}/${Species}" --expected_cell_num "${cell_num}" \
    --cell_calling_method auto --bam "${sampledir}/04.featureCounts/${sample}_nameSorted.bam" --force_cell_num None

# Run mode-specific branch
if [[ "$mode" == "scrna" ]]; then
    echo "scRNA mode: Skipping spatial analysis steps (binSegment, spatial analysis, and spatial report)"
    echo "Pipeline completed successfully at $(date)"
    echo "scRNA analysis results available in: ${sampledir}"
else
    # Run method-specific branch for spatial modes
    if [ "$method" == "image" ]; then
        run_image_branch
    elif [ "$method" == "gene_expr" ]; then
        run_gene_expr_branch
    fi

    ## 08.analysis
    run_command celatlas_spatial rna analysis --outdir "${sampledir}/07.analysis" --sample "${sample}" \
        --thread "${thread}" --genomeDir "${reference_dir}/${Species}" --square_bin_dir "${sampledir}/06.binSegment/square_bin" \
        --pixel-size "${pixelSize}" --bin "${bin}"

    ## Generate comprehensive spatial analysis report
    generate_spatial_report

    echo "Pipeline completed successfully at $(date)"
    echo "Spatial analysis report available at: ${sampledir}/${sample}_spatial_analysis_report.html"
fi
