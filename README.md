# Genome Assembly Tools
This repository contains tools for genome assembly and analysis using different algorithms and approaches. The repository includes:

1. De Bruijn Graph (DBG) genome assembler
2. Overlap-Layout-Consensus (OLC) genome assembler
3. HiFiasm genome assembly pipeline for long reads
4. Assembly evaluation and comparison scripts

## Prerequisites

### For Python Assemblers (DBG and OLC)
- Python 3.6+
- Biopython (`pip install biopython`)
- NetworkX (`pip install networkx`)
- NumPy (`pip install numpy`)
- QUAST (installed and in your PATH)

Install Python dependencies:
```bash
pip install biopython networkx numpy
```

### For Slurm Scripts (HPC Environment)
- Access to an HPC environment with Slurm workload manager
- Required modules: hifiasm, merqury, busco, quast, inspector

## Usage Instructions

### De Bruijn Graph (DBG) Assembler

The DBG assembler creates contigs by building a graph from k-mers extracted from sequencing reads.

```bash
python DeBruijnGraph.py --input <input.fastq> --kmer_size <k> --output <output.fasta> [--reference <reference.fasta>] [--min_coverage <min_cov>] [--export_gfa <graph.gfa>]
```

Example:
```bash
python DeBruijnGraph.py --input reads.fastq --kmer_size 21 --output dbg_assembly.fasta --min_coverage 2
```

Parameters:
- `--input`: Path to the input FASTQ file
- `--kmer_size`: Size of k-mers to use for graph construction
- `--output`: Path to the output FASTA file for contigs
- `--reference`: (Optional) Path to the reference genome for evaluation with QUAST
- `--min_coverage`: (Optional) Minimum k-mer coverage (default: 1)
- `--export_gfa`: (Optional) Export assembly graph in GFA format

### Overlap-Layout-Consensus (OLC) Assembler

The OLC assembler creates contigs by finding overlaps between reads, constructing a graph, and generating consensus sequences.

```bash
python OverlapLayoutConsensus.py --input <input.fastq> --min_overlap <n> --output <output.fasta> [--reference <reference.fasta>] [--min_coverage <min_cov>] [--export_gfa <graph.gfa>]
```

Example:
```bash
python OverlapLayoutConsensus.py --input reads.fastq --min_overlap 50 --output olc_assembly.fasta
```

Parameters:
- `--input`: Path to the input FASTQ file
- `--min_overlap`: Minimum overlap length between reads
- `--output`: Path to the output FASTA file for contigs
- `--reference`: (Optional) Path to the reference genome for evaluation with QUAST
- `--min_coverage`: (Optional) Minimum read coverage for consensus (default: 1)
- `--export_gfa`: (Optional) Export assembly graph in GFA format

### HiFiasm Genome Assembly on HPC (Slurm)

The `hifiasm.slurm` script runs the HiFiasm assembler on a high-performance computing cluster using Slurm.

```bash
sbatch hifiasm.slurm
```

This script:
1. Combines multiple HiFi PacBio reads
2. Decompresses Hi-C reads if needed
3. Runs HiFiasm with Hi-C for phasing
4. Converts the output GFA files to FASTA format

Modify the script variables if needed to point to your input files:
- `COMP_HIFI_1`, `COMP_HIFI_2`, `COMP_HIFI_3`: Compressed PacBio HiFi read files
- `COMP_HIC_1`, `COMP_HIC_2`: Compressed Hi-C read files
- `OUTPUT_DIR`, `TEMP_DIR`: Output and temporary directories

### Lizard Assembly Analysis on HPC (Slurm)

The `LizardAssemblyAnalysis.slurm` script evaluates genome assemblies using multiple tools:

```bash
sbatch LizardAssemblyAnalysis.slurm
```

This script runs:
1. Merqury for k-mer assessment and QV scoring
2. BUSCO for gene completeness assessment
3. QUAST for assembly statistics
4. Inspector for misassembly detection

**Important**: Run `hifiasm.slurm` first to generate the assemblies, and then run `LizardAssemblyAnalysis.slurm` to evaluate them.

## Comparing Multiple Assemblies with QUAST

After generating assemblies with different methods, you can compare them using QUAST:

```bash
# Create directory for QUAST results
mkdir -p ~/cs249_hw2_results/quast_comparison

# Run QUAST with multiple assemblies
python3 /path/to/quast-5.2.0/quast.py \
  -r /path/to/reference/genome.fna \
  -o ~/cs249_hw2_results/quast_comparison \
  -m 100 \
  -l "DBG_Illumina_NoErr,OLC_Illumina_NoErr,DBG_Illumina_Err,OLC_Illumina_Err,DBG_ONT_NoErr,OLC_ONT_NoErr,DBG_ONT_Err,OLC_ONT_Err" \
  /path/to/dbg_hiseq_no_error.fasta \
  /path/to/olc_hiseq_no_error.fasta \
  /path/to/dbg_hiseq_with_error.fasta \
  /path/to/olc_hiseq_with_error.fasta \
  /path/to/dbg_ont_no_error.fasta \
  /path/to/olc_ont_no_error.fasta \
  /path/to/dbg_ont_with_error.fasta \
  /path/to/olc_ont_with_error.fasta
```

For comparing with SPAdes assemblies:

```bash
mkdir -p ~/cs249_hw2_results/quast_spades_comparison

python3 /path/to/quast-5.2.0/quast.py \
  -r /path/to/reference/genome.fna \
  -o ~/cs249_hw2_results/quast_spades_comparison \
  -m 100 \
  -l "DBG_Illumina_NoErr,OLC_Illumina_NoErr,SPAdes_Illumina_NoErr,DBG_Illumina_Err,OLC_Illumina_Err,SPAdes_Illumina_Err,DBG_ONT_NoErr,OLC_ONT_NoErr,DBG_ONT_Err,OLC_ONT_Err,SPAdes_Hybrid_NoErr,SPAdes_Hybrid_Err" \
  /path/to/dbg_hiseq_no_error.fasta \
  /path/to/olc_hiseq_no_error.fasta \
  /path/to/spades_hiseq_no_error/contigs.fasta \
  /path/to/dbg_hiseq_with_error.fasta \
  /path/to/olc_hiseq_with_error.fasta \
  /path/to/spades_hiseq_with_error/contigs.fasta \
  /path/to/dbg_ont_no_error.fasta \
  /path/to/olc_ont_no_error.fasta \
  /path/to/dbg_ont_with_error.fasta \
  /path/to/olc_ont_with_error.fasta \
  /path/to/spades_hybrid_no_error/contigs.fasta \
  /path/to/spades_hybrid_with_error/contigs.fasta
```

## Example WSL Workflow

If working with Windows Subsystem for Linux (WSL), here's a complete workflow:

```bash
# 1. Run the DBG assembler
python DeBruijnGraph.py --input /mnt/c/Users/username/data/reads.fastq --kmer_size 21 --output /mnt/c/Users/username/results/dbg_assembly.fasta

# 2. Run the OLC assembler
python OverlapLayoutConsensus.py --input /mnt/c/Users/username/data/reads.fastq --min_overlap 50 --output /mnt/c/Users/username/results/olc_assembly.fasta

# 3. Compare assemblies with QUAST
mkdir -p ~/cs249_hw2_results/quast_comparison
python3 /path/to/quast-5.2.0/quast.py \
  -r /mnt/c/Users/username/data/reference.fna \
  -o ~/cs249_hw2_results/quast_comparison \
  -m 100 \
  -l "DBG,OLC" \
  /mnt/c/Users/username/results/dbg_assembly.fasta \
  /mnt/c/Users/username/results/olc_assembly.fasta

# 4. Copy results back to Windows if needed
cp -r ~/cs249_hw2_results/quast_comparison /mnt/c/Users/username/results/
```

## Notes

- For the Python assemblers, QUAST must be installed and available in your PATH for the assembly evaluation feature to work.
- When running on large genomes, increase the memory and time allocations in the Slurm scripts or limit the memory on merqury for the analysis.
