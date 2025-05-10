#!/bin/bash -l
#SBATCH --job-name=hifiasm_lizard
#SBATCH --time=48:00:00
#SBATCH --partition=batch
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=256G
#SBATCH --output=hifiasm_lizard.%J.out
#SBATCH --error=hifiasm_lizard.%J.err
#SBATCH --account=cs249

# There are some lines to check if I already have the combined HIFI file since this wasn't the first try and didn't want to combine multiple timestamps
# Similar thing with decompressed files. However, I realized after emailing/meeting with the TAs that I didn't need to decompress.

# Load required modules
module load hifiasm

# Create output directory
OUTPUT_DIR="/ibex/scratch/rajehyl/assembly_results/hifiasm"
TEMP_DIR="/ibex/scratch/rajehyl/temp_decompressed_files"
mkdir -p ${OUTPUT_DIR}

# Define input data paths
COMP_HIFI_1="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver.fastq.gz"
COMP_HIFI_2="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_rev.fastq.gz"
COMP_HIFI_3="/ibex/reference/course/cs249/lizard/input/pacbio/lizard_liver_seq.fastq.gz"
COMP_HIC_1="/ibex/reference/course/cs249/lizard/input/hic/lizard_hic_R1.fastq.gz"
COMP_HIC_2="/ibex/reference/course/cs249/lizard/input/hic/lizard_hic_R2.fastq.gz"

# Define the combined HiFi file path
COMBINED_HIFI="${TEMP_DIR}/hifi/combined_lizard_hifi.fastq"
HIC_1="${TEMP_DIR}/hic/lizard_hic_R1.fastq"
HIC_2="${TEMP_DIR}/hic/lizard_hic_R2.fastq"

# Check if the combined HiFi file already exists
if [ ! -f "${COMBINED_HIFI}" ]; then
    echo "Combined HiFi file not found. Creating it now..."
    mkdir -p ${TEMP_DIR}/hifi

    # Decompress and combine the HiFi files
    srun --ntasks=1 --cpus-per-task=1 sh -c "gzip -dc ${COMP_HIFI_1} > ${COMBINED_HIFI}"
    srun --ntasks=1 --cpus-per-task=1 sh -c "gzip -dc ${COMP_HIFI_2} >> ${COMBINED_HIFI}"
    srun --ntasks=1 --cpus-per-task=1 sh -c "gzip -dc ${COMP_HIFI_3} >> ${COMBINED_HIFI}"

    # Print file size
    echo "Created combined HiFi file:"
    ls -lah ${COMBINED_HIFI}
else
    echo "Using existing combined HiFi file: ${COMBINED_HIFI}"
    ls -lah ${COMBINED_HIFI}
fi

# Check if Hi-C files already exist, if not decompress them
if [ ! -f "${HIC_1}" ]; then
    echo "Hi-C R1 file not found. Decompressing..."
    mkdir -p ${TEMP_DIR}/hic
    srun --ntasks=1 --cpus-per-task=1 gzip -dc ${COMP_HIC_1} > ${HIC_1}
    echo "Created Hi-C R1 file:"
    ls -lah ${HIC_1}
else
    echo "Using existing Hi-C R1 file: ${HIC_1}"
    ls -lah ${HIC_1}
fi

if [ ! -f "${HIC_2}" ]; then
    echo "Hi-C R2 file not found. Decompressing..."
    mkdir -p ${TEMP_DIR}/hic
    srun --ntasks=1 --cpus-per-task=1 gzip -dc ${COMP_HIC_2} > ${HIC_2}
    echo "Created Hi-C R2 file:"
    ls -lah ${HIC_2}
else
    echo "Using existing Hi-C R2 file: ${HIC_2}"
    ls -lah ${HIC_2}
fi

# Run hifiasm with Hi-C data
echo "Running hifiasm with Hi-C data integration..."
OUTPUT_PREFIX="${OUTPUT_DIR}/lizard_hifiasm"

# Run hifiasm with --h1 and --h2 for Hi-C based phasing
srun --ntasks=1 --cpus-per-task=${SLURM_CPUS_PER_TASK} \
    hifiasm \
    -o ${OUTPUT_PREFIX} \
    -t ${SLURM_CPUS_PER_TASK} \
    --h1 ${HIC_1} \
    --h2 ${HIC_2} \
    ${COMBINED_HIFI}

echo "Converting GFA to FASTA format..."
# Convert the output GFA files to FASTA format
if [ -f "${OUTPUT_PREFIX}.hap1.p_ctg.gfa" ]; then
    awk '/^S/{print ">"$2"\n"$3}' ${OUTPUT_PREFIX}.hap1.p_ctg.gfa > ${OUTPUT_PREFIX}.hap1.p_ctg.fa
    echo "Created FASTA for haplotype 1 primary contigs"
fi

if [ -f "${OUTPUT_PREFIX}.hap2.p_ctg.gfa" ]; then
    awk '/^S/{print ">"$2"\n"$3}' ${OUTPUT_PREFIX}.hap2.p_ctg.gfa > ${OUTPUT_PREFIX}.hap2.p_ctg.fa
    echo "Created FASTA for haplotype 2 primary contigs"
fi

if [ -f "${OUTPUT_PREFIX}.p_ctg.gfa" ]; then
    awk '/^S/{print ">"$2"\n"$3}' ${OUTPUT_PREFIX}.p_ctg.gfa > ${OUTPUT_PREFIX}.p_ctg.fa
    echo "Created FASTA for primary contigs"
fi

if [ -f "${OUTPUT_PREFIX}.a_ctg.gfa" ]; then
    awk '/^S/{print ">"$2"\n"$3}' ${OUTPUT_PREFIX}.a_ctg.gfa > ${OUTPUT_PREFIX}.a_ctg.fa
    echo "Created FASTA for alternate contigs"
fi

echo "Hifiasm assembly completed. Results are in ${OUTPUT_PREFIX}.*"
echo "The following output files were created:"
ls -lah ${OUTPUT_PREFIX}*
