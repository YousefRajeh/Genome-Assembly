#!/bin/bash -l
#SBATCH --job-name=eval_assembly
#SBATCH --time=48:00:00
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --output=eval_assembly.%J.out
#SBATCH --error=eval_assembly.%J.err
#SBATCH --account=cs249

# This script evaluates a genome assembly using Merqury, BUSCO, QUAST, and Inspector
# Order of execution: Merqury → BUSCO → QUAST → Inspector
# It takes a long time for Merqury so although it's first here, I actually ran it with QUAST → Inspector → Busco → Merqury

# Function to print status messages with timestamps
log_status() {
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] $1"
}

# Function to check if a file exists and exit if it doesn't
check_file() {
    if [ ! -f "$1" ]; then
        log_status "ERROR: Required file not found: $1"
        exit 1
    else
        log_status "Found file: $1"
    fi
}

# Base directories
ASSEMBLY_DIR="/ibex/scratch/rajehyl/assembly_results/hifiasm"
EVAL_DIR="/ibex/scratch/rajehyl/assembly_results/evaluation_efficient"
# Use /tmp for temporary files (automatically cleaned up when job ends)
TEMP_DIR="/tmp/${USER}_eval_${SLURM_JOB_ID}"
RAW_DATA_DIR="/ibex/reference/course/cs249/lizard/input"

# Create output directories
log_status "Creating output directories..."
mkdir -p ${EVAL_DIR}/{quast,busco,merqury,inspector}
# Create temporary directory in /tmp (will be faster and auto-cleaned)
mkdir -p ${TEMP_DIR}
log_status "Using temporary directory: ${TEMP_DIR}"

# Input files - paths to FASTA assembly files
HAP1_FASTA="${ASSEMBLY_DIR}/lizard_hifiasm.hic.hap1.p_ctg.fa"
HAP2_FASTA="${ASSEMBLY_DIR}/lizard_hifiasm.hic.hap2.p_ctg.fa"

# Verify assembly files exist and are in FASTA format
log_status "Checking for FASTA files..."
for FASTA_FILE in "${HAP1_FASTA}" "${HAP2_FASTA}"; do
    check_file "${FASTA_FILE}"

    # Check if file is in FASTA format (starts with >)
    if [ "$(head -c 1 ${FASTA_FILE})" != ">" ]; then
        log_status "ERROR: File doesn't appear to be in FASTA format: ${FASTA_FILE}"
        log_status "First character should be '>' but is '$(head -c 1 ${FASTA_FILE})'"
        exit 1
    fi
done

# Path to existing combined compressed HiFi file
COMBINED_HIFI_GZ="/ibex/scratch/rajehyl/tmp_hifi_test1/combined_hifi.fastq.gz"
if [ ! -f "${COMBINED_HIFI_GZ}" ]; then
    log_status "Warning: Combined compressed HiFi file not found at: ${COMBINED_HIFI_GZ}"
    log_status "Looking for alternative paths..."

    # Check other possible locations
    for alternative in "/ibex/scratch/rajehyl/tmp_hifi_test1/combined_lizard_hifi.fastq.gz" \
                       "/ibex/scratch/rajehyl/tmp_hifi/combined_hifi.fastq.gz" \
                       "/ibex/scratch/rajehyl/tmp_hifi/combined_lizard_hifi.fastq.gz"; do
        if [ -f "${alternative}" ]; then
            COMBINED_HIFI_GZ="${alternative}"
            log_status "Found combined compressed HiFi file at: ${COMBINED_HIFI_GZ}"
            break
        fi
    done

    # If still not found, use original compressed files as fallback
    if [ ! -f "${COMBINED_HIFI_GZ}" ]; then
        log_status "Could not find existing combined compressed HiFi file."
        log_status "Will use original compressed files instead."
        COMPRESSED_HIFI_1="${RAW_DATA_DIR}/pacbio/lizard_liver.fastq.gz"
        COMPRESSED_HIFI_2="${RAW_DATA_DIR}/pacbio/lizard_liver_rev.fastq.gz"
        COMPRESSED_HIFI_3="${RAW_DATA_DIR}/pacbio/lizard_liver_seq.fastq.gz"

        # Check if the original files exist
        check_file "${COMPRESSED_HIFI_1}"
        check_file "${COMPRESSED_HIFI_2}"
        check_file "${COMPRESSED_HIFI_3}"

        USE_COMBINED=false
    else
        USE_COMBINED=true
    fi
else
    log_status "Found combined compressed HiFi file at: ${COMBINED_HIFI_GZ}"
    USE_COMBINED=true
fi

log_status "====================================================================="
log_status "Starting memory-efficient assembly evaluation for lizard genome"
log_status "Haplotype 1 assembly: ${HAP1_FASTA}"
log_status "Haplotype 2 assembly: ${HAP2_FASTA}"
if [ "${USE_COMBINED}" = true ]; then
    log_status "Using existing combined compressed HiFi file: ${COMBINED_HIFI_GZ}"
else
    log_status "Using original compressed files"
fi
log_status "====================================================================="

# Install matplotlib for QUAST plotting
log_status "Installing matplotlib to enable plotting in QUAST..."
pip install matplotlib --user

# Decompress the reads to /tmp for faster processing
log_status "Decompressing reads to temporary directory for faster processing..."
TEMP_READS="${TEMP_DIR}/combined_reads.fastq"

# Use either the combined file or the individual files
if [ "${USE_COMBINED}" = true ]; then
    log_status "Decompressing combined file (this may take some time)..."
    zcat ${COMBINED_HIFI_GZ} > ${TEMP_READS}
else
    log_status "Decompressing and combining individual files (this may take some time)..."
    zcat ${COMPRESSED_HIFI_1} ${COMPRESSED_HIFI_2} ${COMPRESSED_HIFI_3} > ${TEMP_READS}
fi

# Check if decompression was successful
if [ ! -f "${TEMP_READS}" ]; then
    log_status "ERROR: Failed to create decompressed reads file."
    exit 1
fi

log_status "Decompression complete. Temporary reads file size: $(du -h ${TEMP_READS} | cut -f1)"

# 1. (FIRST) Merqury: K-mer distribution and QV score - WITH MEMORY LIMIT FIXED
log_status "STEP 1: Running Merqury with memory limit..."
module load merqury

# Change to the merqury output directory
cd ${EVAL_DIR}/merqury

# Clean up any existing merqury directory to avoid errors
if [ -d "${EVAL_DIR}/merqury/read_db.meryl" ]; then
    log_status "Removing existing meryl database to start fresh..."
    rm -rf ${EVAL_DIR}/merqury/read_db.meryl
fi

# Build k-mer database with memory limit (170GB) - FIX FOR MEMORY ISSUE
log_status "Building k-mer database from reads with 170GB memory limit..."
log_status "Command: meryl count k=21 memory=170 ${TEMP_READS} output ${EVAL_DIR}/merqury/read_db.meryl"
meryl count k=21 memory=170 ${TEMP_READS} output ${EVAL_DIR}/merqury/read_db.meryl

# Check if the k-mer database was created successfully
if [ ! -d "${EVAL_DIR}/merqury/read_db.meryl" ]; then
    log_status "ERROR: Failed to create k-mer database."
    exit 1
fi

log_status "K-mer database created successfully."

# Run Merqury on both haplotypes
log_status "Running Merqury on both haplotypes..."
log_status "Command: $MERQURY/merqury.sh ${EVAL_DIR}/merqury/read_db.meryl ${HAP1_FASTA} ${HAP2_FASTA} ${EVAL_DIR}/merqury/lizard_assembly"
$MERQURY/merqury.sh \
    ${EVAL_DIR}/merqury/read_db.meryl \
    ${HAP1_FASTA} ${HAP2_FASTA} \
    ${EVAL_DIR}/merqury/lizard_assembly

# Check if the QV file was created
if [ -f "${EVAL_DIR}/merqury/lizard_assembly.qv" ]; then
    log_status "Merqury completed successfully. QV file created."
    log_status "QV scores: $(cat ${EVAL_DIR}/merqury/lizard_assembly.qv)"
else
    log_status "WARNING: Merqury QV file not found. Merqury may have failed."
fi

log_status "Merqury completed. Results in ${EVAL_DIR}/merqury"

# 2. BUSCO: Gene completeness with auto-lineage - FIX FOR LINEAGE ISSUE
log_status "STEP 2: Running BUSCO with auto-lineage for gene completeness assessment..."
module load busco

# Force flag to overwrite existing results
FORCE_FLAG="-f"

# Clean up any existing busco directories to avoid errors
for HAP_NUM in 1 2; do
    if [ -d "${EVAL_DIR}/busco/busco_hap${HAP_NUM}" ]; then
        log_status "Removing existing BUSCO directory for haplotype ${HAP_NUM} to start fresh..."
        rm -rf ${EVAL_DIR}/busco/busco_hap${HAP_NUM}
    fi
done

# Create a directory for BUSCO downloads
BUSCO_DOWNLOADS="${EVAL_DIR}/busco/downloads"
mkdir -p ${BUSCO_DOWNLOADS}

# Run BUSCO on each haplotype separately WITH AUTO-LINEAGE
for HAP_NUM in 1 2; do
    HAP_FASTA="${ASSEMBLY_DIR}/lizard_hifiasm.hic.hap${HAP_NUM}.p_ctg.fa"

    # Verify the file exists one more time
    check_file "${HAP_FASTA}"

    log_status "Running BUSCO on haplotype ${HAP_NUM} with auto-lineage..."
    log_status "Command: busco -i ${HAP_FASTA} -o busco_hap${HAP_NUM} -m genome --auto-lineage -c ${SLURM_CPUS_PER_TASK} ${FORCE_FLAG} --out_path ${EVAL_DIR}/busco --download_path ${BUSCO_DOWNLOADS}"

    busco \
        -i ${HAP_FASTA} \
        -o busco_hap${HAP_NUM} \
        -m genome \
        --auto-lineage \
        -c ${SLURM_CPUS_PER_TASK} \
        ${FORCE_FLAG} \
        --out_path ${EVAL_DIR}/busco \
        --download_path ${BUSCO_DOWNLOADS}

    # Check if BUSCO completed successfully
    BUSCO_SUMMARY=$(find ${EVAL_DIR}/busco/busco_hap${HAP_NUM} -name "short_summary.*.txt" | head -n 1)
    if [ -n "${BUSCO_SUMMARY}" ] && [ -f "${BUSCO_SUMMARY}" ]; then
        log_status "BUSCO completed successfully for haplotype ${HAP_NUM}."
        log_status "Summary file: ${BUSCO_SUMMARY}"
    else
        log_status "WARNING: BUSCO summary file not found for haplotype ${HAP_NUM}."
    fi
done

log_status "BUSCO completed. Results in ${EVAL_DIR}/busco"

# 3. QUAST: Basic assembly metrics - Now with plots enabled
log_status "STEP 3: Running QUAST for basic assembly metrics..."
module load quast

# Run QUAST with plots enabled (matplotlib now installed)
log_status "Command: quast.py -o ${EVAL_DIR}/quast -t ${SLURM_CPUS_PER_TASK} -e --large --min-contig 1000 --circos --space-efficient --labels Haplotype1,Haplotype2 ${HAP1_FASTA} ${HAP2_FASTA}"

quast.py \
    -o ${EVAL_DIR}/quast \
    -t ${SLURM_CPUS_PER_TASK} \
    -e \
    --large \
    --min-contig 1000 \
    --circos \
    --space-efficient \
    --labels "Haplotype1,Haplotype2" \
    ${HAP1_FASTA} ${HAP2_FASTA}

# Check if QUAST completed successfully
if [ -f "${EVAL_DIR}/quast/report.txt" ]; then
    log_status "QUAST completed successfully."
    log_status "Summary report generated at: ${EVAL_DIR}/quast/report.txt"
else
    log_status "WARNING: QUAST report file not found."
fi

log_status "QUAST completed. Results in ${EVAL_DIR}/quast"

# 4. Inspector: Process a subset for misassembly identification
log_status "STEP 4: Running Inspector on selected longest contigs only..."
module load inspector

# Run Inspector on a subset of contigs to save space
for HAP_NUM in 1 2; do
    HAP_FASTA="${ASSEMBLY_DIR}/lizard_hifiasm.hic.hap${HAP_NUM}.p_ctg.fa"
    OUT_DIR="${EVAL_DIR}/inspector/hap${HAP_NUM}"
    mkdir -p ${OUT_DIR}

    # Create a subset of the largest contigs (first 20) to analyze
    log_status "Creating subset of largest contigs for Haplotype ${HAP_NUM}..."

    # Extract contig names and sizes, sort by size, and keep top 20
    grep "^>" ${HAP_FASTA} | \
    awk '{header=$0; getline; print length($0)"\t"header}' ${HAP_FASTA} | \
    sort -nr | head -n 20 | \
    awk '{print $2}' | sed 's/>//' > ${OUT_DIR}/largest_contigs.txt

    # Create FASTA with only these contigs
    SUBSET_FASTA="${TEMP_DIR}/largest_contigs_hap${HAP_NUM}.fa"
    touch ${SUBSET_FASTA}

    while read CONTIG; do
        awk -v contig="$CONTIG" 'BEGIN {print_flag=0}
        $0 ~ "^>"contig"( |$)" {print_flag=1; print; next}
        $0 ~ "^>" {print_flag=0}
        print_flag {print}' ${HAP_FASTA} >> ${SUBSET_FASTA}
    done < ${OUT_DIR}/largest_contigs.txt

    # Save a copy of the subset in the output directory too
    cp ${SUBSET_FASTA} ${OUT_DIR}/

    log_status "Running Inspector on haplotype ${HAP_NUM}..."
    log_status "Command: inspector.py -c ${SUBSET_FASTA} -r ${TEMP_READS} -d hifi -o ${OUT_DIR} -t ${SLURM_CPUS_PER_TASK} --min_contig_length 50000 --min_contig_length_assemblyerror 100000 --skip_base_error --noplot"

    # Run Inspector with the decompressed reads file
    inspector.py \
        -c ${SUBSET_FASTA} \
        -r ${TEMP_READS} \
        -d hifi \
        -o ${OUT_DIR} \
        -t ${SLURM_CPUS_PER_TASK} \
        --min_contig_length 50000 \
        --min_contig_length_assemblyerror 100000 \
        --skip_base_error \
        --noplot

    # Check if Inspector completed successfully
    if [ -f "${OUT_DIR}/summary.stats" ]; then
        log_status "Inspector completed successfully for haplotype ${HAP_NUM}."
        log_status "Summary stats generated at: ${OUT_DIR}/summary.stats"
    else
        log_status "WARNING: Inspector summary stats not found for haplotype ${HAP_NUM}."
    fi
done

log_status "Inspector completed on subset of longest contigs. Results in ${EVAL_DIR}/inspector"

# Create summary report
SUMMARY_FILE="${EVAL_DIR}/assembly_evaluation_summary.txt"

log_status "Generating summary report..."

{
    echo "======================================================================="
    echo "ASSEMBLY EVALUATION SUMMARY FOR LIZARD GENOME (Scincus mitranus)"
    echo "Date: $(date)"
    echo "======================================================================="
    echo ""

    # Merqury summary (first since it was run first)
    echo "MERQURY QV SCORES AND K-MER COMPLETENESS:"
    echo "-----------------------------------------"
    if [ -f "${EVAL_DIR}/merqury/lizard_assembly.qv" ]; then
        echo "QV Scores:"
        cat ${EVAL_DIR}/merqury/lizard_assembly.qv
    else
        echo "Merqury QV report not found"
    fi

    if [ -f "${EVAL_DIR}/merqury/lizard_assembly.completeness.stats" ]; then
        echo "K-mer Completeness Stats:"
        cat ${EVAL_DIR}/merqury/lizard_assembly.completeness.stats
    else
        echo "Merqury completeness stats not found"
    fi
    echo ""

    # BUSCO summary
    echo "BUSCO GENE COMPLETENESS (Auto-lineage):"
    echo "----------------------------------------"
    echo "Haplotype 1:"
    BUSCO_HAP1=$(find ${EVAL_DIR}/busco/busco_hap1 -name "short_summary.*.txt" | head -n 1)
    if [ -n "${BUSCO_HAP1}" ] && [ -f "${BUSCO_HAP1}" ]; then
        cat "${BUSCO_HAP1}"
    else
        echo "BUSCO Hap1 report not found"
    fi
    echo ""

    echo "Haplotype 2:"
    BUSCO_HAP2=$(find ${EVAL_DIR}/busco/busco_hap2 -name "short_summary.*.txt" | head -n 1)
    if [ -n "${BUSCO_HAP2}" ] && [ -f "${BUSCO_HAP2}" ]; then
        cat "${BUSCO_HAP2}"
    else
        echo "BUSCO Hap2 report not found"
    fi
    echo ""

    # QUAST summary
    echo "QUAST BASIC METRICS:"
    echo "--------------------"
    if [ -f "${EVAL_DIR}/quast/report.txt" ]; then
        cat ${EVAL_DIR}/quast/report.txt
    else
        echo "QUAST report not found"
    fi
    echo ""

    # Inspector summary
    echo "INSPECTOR MIS-ASSEMBLY DETECTION (Top 20 largest contigs only, as it was taking too much space and I ran out multiple times on IBEX):"
    echo "------------------------------------------------------------"
    for HAP_NUM in 1 2; do
        echo "Haplotype ${HAP_NUM}:"

        # QV summary
        if [ -f "${EVAL_DIR}/inspector/hap${HAP_NUM}/summary.stats" ]; then
            echo "QV Summary:"
            cat ${EVAL_DIR}/inspector/hap${HAP_NUM}/summary.stats
        else
            echo "Inspector QV summary not found for Hap${HAP_NUM}"
        fi

        # Structural errors
        if [ -f "${EVAL_DIR}/inspector/hap${HAP_NUM}/structural_error.bed" ]; then
            STRUCT_ERRORS=$(wc -l < ${EVAL_DIR}/inspector/hap${HAP_NUM}/structural_error.bed)
            echo "Number of structural errors detected: ${STRUCT_ERRORS}"

            if [ "$STRUCT_ERRORS" -gt 0 ] && [ "$STRUCT_ERRORS" -lt 20 ]; then
                echo "Structural errors:"
                cat ${EVAL_DIR}/inspector/hap${HAP_NUM}/structural_error.bed
            elif [ "$STRUCT_ERRORS" -ge 20 ]; then
                echo "First 20 structural errors (of ${STRUCT_ERRORS} total):"
                head -n 20 ${EVAL_DIR}/inspector/hap${HAP_NUM}/structural_error.bed
            fi
        else
            echo "No structural errors file found for Hap${HAP_NUM}"
        fi
        echo ""
    done

    echo "======================================================================="
    echo "End of summary"
    echo "For detailed results, see the output directories for each tool."
    echo "======================================================================="

} > ${SUMMARY_FILE}

# No need to manually clean up /tmp - it will be automatically cleaned when the job ends
log_status "Temporary files in ${TEMP_DIR} will be automatically cleaned up when the job ends"

log_status "Evaluation complete! Summary report available at ${SUMMARY_FILE}"
log_status "====================================================================="
log_status "Results directories:"
log_status "Merqury:   ${EVAL_DIR}/merqury"
log_status "BUSCO:     ${EVAL_DIR}/busco"
log_status "QUAST:     ${EVAL_DIR}/quast"
log_status "Inspector: ${EVAL_DIR}/inspector (limited to largest contigs)"
log_status "====================================================================="