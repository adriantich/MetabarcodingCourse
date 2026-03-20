# Curso de Metabarcoding - Mallorca 2025

## Course Description

Welcome to the Metabarcoding Course repository! This course provides hands-on training in metabarcoding bioinformatics, focusing on the analysis of COI mitochondrial gene amplicon sequencing data for metazoan community characterization.

Metabarcoding is a powerful molecular technique that combines DNA barcoding with high-throughput sequencing to identify and quantify organisms in complex environmental samples. This course covers the complete bioinformatics workflow from raw sequencing reads to taxonomic identification and diversity analysis.

## Repository Structure

```
CursoMetabarcodingMallorca2025/
├── data/                          # Raw FASTQ files (paired-end)
├── scripts/                       # Analysis scripts
│   └── pipeline_explained.sh      # Main pipeline execution script
├── tools/                         # Additional tools
├── README.md                      # This file
└── LICENSE                        # GPL-3.0 License
```

## Software Requirements

### Core Tools

The metabarcoding pipeline requires the following software:

1. **FastQC** (≥ 0.11.8) - Quality control of raw sequencing data
2. **Cutadapt** (≥ 5.2) - Adapter trimming and quality filtering
3. **VSEARCH** (≥ 2.30.1) - Sequence analysis, merging, and clustering
4. **dnoise** (≥ 1.4.2) - Denoising of amplicon sequences
5. **SWARM** (≥ 3.1.6) - Clustering algorithm for OTU generation
6. **BLAST+** (≥ 2.17.0) - Taxonomic assignment
7. **R** (≥ 4.4.3) - Statistical analysis and MJOLNIR3 package
8. **mumu** - Post-clustering curation tool
9. **mkLTG** - Local taxonomy generator

### Required R Packages

- **MJOLNIR3** - Main metabarcoding pipeline framework
- **Biostrings** (≥ 2.74.0) - DNA sequence manipulation
- **Rcpp** (≥ 1.1.0) - R/C++ interface
- **dplyr** (≥ 1.1.4) - Data manipulation
- **tidyr** (≥ 1.3.1) - Data tidying
- **stringr** (≥ 1.6.0) - String operations

### Taxonomy Reference Database

A taxonomy reference database is required for taxonomic assignment. Download from the provided Google Drive link.

## Installation Instructions

### Option 1: Using Conda/Mamba (Recommended)

Conda provides an easy way to install all required tools in an isolated environment.

#### Step 1: Install Miniconda or Mambaforge

**Download and install Miniconda:**
```bash
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-Linux-x86_64.sh
```

#### Step 2: Create a Conda Environment

```bash

# clone the repository
git clone https://github.com/adriantich/CursoMetabarcodingMallorca2025.git
cd CursoMetabarcodingMallorca2025

# Create environment with all required tools
# conda can also be used but mamba is prefered
mamba create -n metabarcoding -c bioconda -c conda-forge python=3.11.14 \
    fastqc=0.11.8 \
    cutadapt=5.2 \
    vsearch=2.30.1 \
    dnoise=1.4.2 \
    swarm=3.1.6 \
    r-base=4.4.3 \
    bioconductor-biostrings=2.74.0 \
    r-rcpp=1.1.0 \
    r-dplyr=1.1.4 \
    r-tidyr=1.3.1 \
    r-stringr=1.6.0 \
    cxx-compiler=1.0.0

# Activate the environment
mamba activate metabarcoding

mkdir -p SOFT
cd SOFT

wget https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-x64-linux.tar.gz
tar zxvpf ncbi-blast-2.17.0+-x64-linux.tar.gz

cp -r ncbi-blast-2.17.0+/bin/* $CONDA_PREFIX/bin/.

git clone https://github.com/meglecz/mkLTG.git

git clone https://github.com/frederic-mahe/mumu.git
cd mumu
make ; make check ; make install prefix=$CONDA_PREFIX
cd ..
# warning: if g++ not detected a broken symlink problem could be the cause
# remake the link by 
# cd $CONDA_PREFIX/bin && rm g++ && ln -s x86_64-conda-linux-gnu-g++ g++

git clone https://github.com/adriantich/MJOLNIR3.git
cd MJOLNIR3
git checkout obitools2vsearch
cd ..
R CMD INSTALL MJOLNIR3

cd ../tools
make ; make install PREFIX=$CONDA_PREFIX
cd ..

```

### Download taxonomy reference database

#### Manually
https://drive.google.com/drive/folders/1-LBlUKFA-r5g6GI0sTo-7t92ml3pAS0S?usp=sharing

```bash
mamba activate metabarcoding
mamba install conda-forge::gdown
gdown --folder https://drive.google.com/drive/folders/1-LBlUKFA-r5g6GI0sTo-7t92ml3pAS0S
```
### Verification

Test that all tools are properly installed:

```bash
# Check FastQC
fastqc --version
# Expected output: FastQC v0.11.9 (or higher)

# Check Cutadapt
cutadapt --version
# Expected output: 3.5 (or higher)

# Check VSEARCH
vsearch --version
# Expected output: vsearch v2.22.1 (or higher)

# Check BLAST+
blastn -version
# Expected output: blastn: 2.17.0+ (or higher)

# Check dnoise
dnoise --version
# Expected output: dnoise 1.4.2 (or higher)

# Check swarm
swarm --version
# Expected output: Swarm 3.1.6 (or higher)

# Check R
R --version
# Expected output: R version 4.4.3 (or higher)

# Check mumu
mumu --version
# Expected output: mumu version information

# Check MJOLNIR3 in R
R -e "library(mjolnir)"
# Expected output: No errors, package loaded successfully


```

## Quick Start Guide

### 1. Activate Your Environment (if using Conda)

```bash
mamba activate metabarcoding
```

### 3. Run the Pipeline for a single sample as example

```bash
cd scripts
./pipeline_explained.sh
```

The pipeline will:
- Perform quality control on raw sequences
- Demultiplex the samples
- Merge paired-end reads
- Filter and trim reads based on quality
- Dereplicate and remove chimeras
- Denoise to obtain ESV and cluster them into OTUs
- Taxonomic assignment
- Post-clustering filter


## Dataset Information

The `data/` directory contains small demonstration datasets:
- **Sample format:** Paired-end FASTQ files (Illumina)
- **Target region:** 16S rRNA V4 hypervariable region
- **Sequencing platform:** Illumina MiSeq (2x250 bp)
- **Purpose:** Educational demonstration only

For working with your own data, replace the files in `data/raw_sequences/` with your samples following the same naming convention: `samplename_R1.fastq` and `samplename_R2.fastq`.

## Common Issues

**Issue 1: "Command not found" errors**
- Solution: Ensure your conda environment is activated or tools are in your PATH
- Verify installation with `which toolname`

**Issue 2: Pipeline fails while running the software**
- Solution: Check that input files are in valid format
- Verify sufficient disk space is available
- Review software log files for specific errors

**Issue 3: Results do not match the expected results**
- Solution: Check the number of reads at each step to detect were the problem is
- Review the metadata information
- Verify correct parameters adapted to your amplicon

For additional help, please contact course instructors.


## Citation

If you use these materials in your research or teaching, please cite:

```
Metabarcoding Course Materials - Mallorca 2025
https://github.com/adriantich/CursoMetabarcodingMallorca2025
```

## License

This project is licensed under the GNU General Public License v3.0 - see the [LICENSE](LICENSE) file for details.

---

**Last updated:** November 2025  
**Course website:** [25-27 of November 2025]