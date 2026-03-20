#!/bin/bash

################################################################################
# Metabarcoding Pipeline - Mallorca 2025 Course
# This script has the commands to run the complete metabarcoding analysis pipeline
# copy and paste them into your terminal
################################################################################

# We first define the variables that will be used throughout the pipeline

R1_FILES="ULO1_R1.fastq.gz"
LIB_PREFIX="ULO1"
EXPERIMENT="ULOY"
CORES=5

cd data
mkdir -p quality_control

# Step 1: Quality Check
echo "Step 1: Quality Check Analysis"
echo "Running FastQC on raw sequences..."

fastqc ULO1_R1.fastq.gz -o quality_control/ --quiet 2>/dev/null || true
echo "  ✓ Quality check completed"

# open the html reports in your browser to visualize them

################################################################################
################################################################################
# Step 2: Demultiplexing and Primer Removal
echo ""
echo "Step 2: Demultiplexing and Primer Removal"
echo "Running demultiplexing and primer removal using cutadapt..." 
echo "We are going to extract que first sample from the R1 files:"
head -n 2 ULOY_metadata.tsv
# "mjolnir_agnomens"	"original_samples"	"date"	        "group"	"position"	"BLANK"
# "ULO1_sample_006"	    "10_25_M1_A_C"	    "2017_10_25"	"cage"	"cage_M1"
echo ""
echo "ngsfile information:"
head -n 1 ngsfilter_ULO1.tsv
# ULO1	ULO1_sample_006	ggatgatc:ggatgatc	GGWACWRGWTGRACWNTNTAYCCYCC	TANACYTCNGGRTGNCCRAARAAYCA
echo ""

# variables for the demultiplexing step
# For one sample
R1_FILES="ULO1_R1.fastq.gz"
R2_FILES="ULO1_R2.fastq.gz"
LIB_PREFIX="ULO1"
EXPERIMENT="ULOY"
CORES=5
ORIGINAL_SAMPLE_NAME="10_25_M1_A_C"
SAMPLE_NAME="ULO1_sample_006"
FWD_tag="ggatgatc"
REV_tag="ggatgatc"
FWD_PRIMER="GGWACWRGWTGRACWNTNTAYCCYCC"
REV_PRIMER="TANACYTCNGGRTGNCCRAARAAYCA"

echo -e ">${SAMPLE_NAME}_sametag\n${FWD_tag}" >fwd_tag.fasta
echo -e ">${SAMPLE_NAME}_sametag\n${REV_tag}" >rev_tag.fasta

# First we demultiplex the sample tags using cutadapt
cutadapt -e 0 -O 8 --no-indels -j $CORES --action='none' --max-n=0.5 \
    --pair-adapters -g file:fwd_tag.fasta -G file:rev_tag.fasta \
    -o {name}_R1.fastq -p {name}_R2.fastq \
    $R1_FILES $R2_FILES
    # allow 0 errors
    # min 8 overlap required
    # no indels allowed
    # number $CORES of cores allowed
    # don't remove the tag taken
    # save those reads that have not been assigned to the sample
    # save those reads that have not been assigned to the sample
    # I allow a max of half of the read being N. this will be solved in following steps

echo -e ">${SAMPLE_NAME}_fwd_rev\n${FWD_PRIMER}\n>${SAMPLE_NAME}_rev_fwd\n${REV_PRIMER}" >R1_primers.fasta
echo -e ">${SAMPLE_NAME}_fwd_rev\n${REV_PRIMER}\n>${SAMPLE_NAME}_rev_fwd\n${FWD_PRIMER}" >R2_primers.fasta

# Then we remove the primers from the demultiplexed files
cutadapt -e 0.1 --no-indels -j $CORES --action='trim' --max-n=0.5 \
    --pair-adapters -g file:R1_primers.fasta -G file:R2_primers.fasta \
    -o {name}_R1.fastq -p {name}_R2.fastq \
    "${SAMPLE_NAME}_sametag_R1.fastq" "${SAMPLE_NAME}_sametag_R2.fastq"
    # allow 10% errors
    # no indels allowed
    # number $CORES of cores allowed
    # remove the primers

cat ${SAMPLE_NAME}_fwd_rev_R1.fastq \
    ${SAMPLE_NAME}_rev_fwd_R2.fastq > ${ORIGINAL_SAMPLE_NAME}_R1.fastq

cat ${SAMPLE_NAME}_fwd_rev_R2.fastq \
    ${SAMPLE_NAME}_rev_fwd_R1.fastq > ${ORIGINAL_SAMPLE_NAME}_R2.fastq

rm ${SAMPLE_NAME}_fwd_rev_R*.fastq ${SAMPLE_NAME}_rev_fwd_R*.fastq \
    ${SAMPLE_NAME}_sametag_R*.fastq unknown*_R*.fastq


################################################################################
# step 2 for all samples
# cutadapt -e 0 -O 8 --no-indels -j $CORES --action='none' --max-n=0.5 \
#     --pair-adapters -g file:tags_all.fasta -G file:tags_all.fasta \
#     -o {name}_R1.fastq -p {name}_R2.fastq \
#     $R1_FILES $R2_FILES

# for R1_FILE in *sample*_R1.fastq; do
#     SAMPLE_NAME=$(basename "$R1_FILE" _R1.fastq)
#     ORIGINAL_SAMPLE_NAME=$(grep "$SAMPLE_NAME" ULOY_metadata.tsv | cut -f2)
#     R2_FILE="${SAMPLE_NAME}_R2.fastq"
#     echo -e ">${SAMPLE_NAME}_fwd_rev\n${FWD_PRIMER}\n>${SAMPLE_NAME}_rev_fwd\n${REV_PRIMER}" >R1_primers.fasta
#     echo -e ">${SAMPLE_NAME}_fwd_rev\n${REV_PRIMER}\n>${SAMPLE_NAME}_rev_fwd\n${FWD_PRIMER}" >R2_primers.fasta

#     cutadapt -e 0.1 --no-indels -j $CORES --action='trim' --max-n=0.5 \
#         --pair-adapters -g file:R1_primers.fasta -G file:R2_primers.fasta \
#         -o {name}_R1.fastq -p {name}_R2.fastq \
#         "${SAMPLE_NAME}_sametag_R1.fastq" "${SAMPLE_NAME}_sametag_R2.fastq"
#     cat ${SAMPLE_NAME}_fwd_rev_R1.fastq ${SAMPLE_NAME}_rev_fwd_R2.fastq > ${ORIGINAL_SAMPLE_NAME}_R1.fastq
#     cat ${SAMPLE_NAME}_fwd_rev_R2.fastq ${SAMPLE_NAME}_rev_fwd_R1.fastq > ${ORIGINAL_SAMPLE_NAME}_R2.fastq
#     rm ${SAMPLE_NAME}_fwd_rev_R*.fastq ${SAMPLE_NAME}_rev_fwd_R*.fastq ${SAMPLE_NAME}_sametag_R*.fastq
# done

# rm unknown*_R*.fastq

################################################################################
################################################################################
# Step 3: Merging Paired-End Reads
echo "Step 3: Merging Paired-End Reads" 
echo "Merging forward and reverse reads..."


LIB_PREFIX="ULO1"
EXPERIMENT="ULOY"
CORES=5
ORIGINAL_SAMPLE_NAME="10_25_M1_A_C"
SAMPLE_NAME="ULO1_sample_006"


vsearch --fastq_mergepairs "${ORIGINAL_SAMPLE_NAME}_R1.fastq" \
    --reverse "${ORIGINAL_SAMPLE_NAME}_R2.fastq" \
    --fastqout "${EXPERIMENT}_${SAMPLE_NAME}_aligned.fastq" \
    --fastq_minovlen 40 \
    --fastq_maxdiffs 0 \
    --fastq_maxee 0.5 \
    --fastq_minmergelen 299 \
    --fastq_maxmergelen 320 \
    --fastq_maxns 0
    # merge pairs of sequences
    # min overlap of 40 bp
    # no differences allowed in the overlap
    # max expected errors of 0.5
    # min merged length of 299 bp
    # max merged length of 320 bp
    # max number of Ns allowed is 0

vsearch --fastx_uniques "${EXPERIMENT}_${SAMPLE_NAME}_aligned.fastq" \
    --sizeout \
    --fastaout "${EXPERIMENT}_${SAMPLE_NAME}_uniq.fasta"
    # dereplicate sequences
    # output in fasta format with size annotation

################################################################################
# step 3 for all samples


# EXPERIMENT="ULOY"

# for R1_FILE in *_R1.fastq; do
#     ORIGINAL_SAMPLE_NAME=$(basename "$R1_FILE" _R1.fastq)
#     SAMPLE_NAME=$(grep "$ORIGINAL_SAMPLE_NAME" ULOY_metadata.tsv | cut -f1)
#     R2_FILE="${ORIGINAL_SAMPLE_NAME}_R2.fastq"
    
#     vsearch --fastq_mergepairs "$R1_FILE" \
#         --reverse "$R2_FILE" \
#         --fastqout "${EXPERIMENT}_${SAMPLE_NAME}_aligned.fastq" \
#         --fastq_minovlen 40 \
#         --fastq_maxdiffs 0 \
#         --fastq_maxee 0.5 \
#         --fastq_minmergelen 299 \
#         --fastq_maxmergelen 320 \
#         --fastq_maxns 0
    
#     vsearch --fastx_uniques "${EXPERIMENT}_${SAMPLE_NAME}_aligned.fastq" \
#         --sizeout \
#         --fastaout "${EXPERIMENT}_${SAMPLE_NAME}_uniq.fasta"

# done

################################################################################
################################################################################
# Step 4: Quimera Removal


LIB_PREFIX="ULO1"
EXPERIMENT="ULOY"
CORES=5
ORIGINAL_SAMPLE_NAME="10_25_M1_A_C"
SAMPLE_NAME="ULO1_sample_006"

vsearch --uchime_denovo "${EXPERIMENT}_${SAMPLE_NAME}_uniq.fasta" \
    --nonchimeras "${EXPERIMENT}_${SAMPLE_NAME}_nochimera.fasta" \
    --chimeras "${EXPERIMENT}_${SAMPLE_NAME}_chimera.fasta" \
    --uchimeout "${EXPERIMENT}_${SAMPLE_NAME}_uchime.log" \
    --threads $CORES \
    --sizeout \
    --minh 0.9
    # de novo chimera detection
    # output non-chimeric sequences
    # output chimeric sequences
    # uchime log file
    # number $CORES of cores allowed
    # size annotation kept
    # minimum score of 0.9 to be considered a chimera
################################################################################
# step 4 for all samples


# for INPUT_FILE in *_uniq.fasta; do
#     SAMPLE_NAME=$(basename "$INPUT_FILE" _uniq.fasta)
#     vsearch --uchime_denovo "$INPUT_FILE" \
#         --nonchimeras "${SAMPLE_NAME}_nochimera.fasta" \
#         --chimeras "${SAMPLE_NAME}_chimera.fasta" \
#         --uchimeout "${SAMPLE_NAME}_uchime.log" \
#         --threads $CORES \
#         --sizeout \
#         --minh 0.9
# done

################################################################################
################################################################################
# Step 5: Obtain ESV and OTU Tables


LIB_PREFIX="ULO1"
EXPERIMENT="ULOY"
CORES=5
ORIGINAL_SAMPLE_NAME="10_25_M1_A_C"
SAMPLE_NAME="ULO1_sample_006"

dnoise --fasta_input "${EXPERIMENT}_${SAMPLE_NAME}_nochimera.fasta" \
    --fasta_output "${EXPERIMENT}_${SAMPLE_NAME}" \
    -c $CORES \
    -a 4 \
    -y -e 0.47,0.23,1.02 \
    -m 313 \
    -r 1
    # denoise sequences to obtain ESVs
    # output ESV sequences
    # number $CORES of cores allowed
    # use alpha of 4 (see Edgar 2016)
    # use entropy correction with these values
    # modal length of 313 bp
    # minimum size of 1

# annotate ESVs with sample name
sed "s/;$/;sample=${EXPERIMENT}_${SAMPLE_NAME}/g" \
    "${EXPERIMENT}_${SAMPLE_NAME}_Adcorr_denoised_ratio_d.fasta" > \
    "${EXPERIMENT}_${SAMPLE_NAME}_annotated.fasta"

## to install rename_fasta go to tools:
# make 
# make install PREFIX=$CONDA_PREFIX

rename_fasta "${EXPERIMENT}_${SAMPLE_NAME}_annotated.fasta" \
    "${EXPERIMENT}_${SAMPLE_NAME}_renamed.fasta" \
    "${EXPERIMENT}"
# simplify headers names

# remove sample name from headers
sed -i "s/;sample=${EXPERIMENT}_${SAMPLE_NAME}//g" \
    "${EXPERIMENT}_${SAMPLE_NAME}_renamed.fasta"

# create ESV table
vsearch --search_exact "${EXPERIMENT}_${SAMPLE_NAME}_annotated.fasta" \
    --db "${EXPERIMENT}_${SAMPLE_NAME}_renamed.fasta" \
    --otutabout "${EXPERIMENT}_${SAMPLE_NAME}_esv_table.csv" \
    --threads $CORES
    # create ESV table
    # number $CORES of cores allowed


swarm -d 13 \
    -z \
    -t $CORES \
    -o "${EXPERIMENT}_${SAMPLE_NAME}_otu_table.txt" \
    -s "${EXPERIMENT}_${SAMPLE_NAME}_swarm_stats.txt" \
    -w "${EXPERIMENT}_${SAMPLE_NAME}_otus.fasta" \
    "${EXPERIMENT}_${SAMPLE_NAME}_renamed.fasta"

sed -i 's/;.*//g' "${EXPERIMENT}_${SAMPLE_NAME}_renamed.fasta"

seq2tab "${EXPERIMENT}_${SAMPLE_NAME}_esv_table.csv" \
    "${EXPERIMENT}_${SAMPLE_NAME}_renamed.fasta" \
    "ID"
    # create seq2tab file

assign_motu \
    "${EXPERIMENT}_${SAMPLE_NAME}_otu_table.txt" \
    "${EXPERIMENT}_${SAMPLE_NAME}_esv_table.csv"

Rscript -e '
library(dplyr)
esv_table <- read.csv("'${EXPERIMENT}'_'${SAMPLE_NAME}'_esv_table.csv", header = TRUE, sep = "\t")
esv_table_sum <- esv_table %>%
    select(MOTU, starts_with("'${EXPERIMENT}'_'${SAMPLE_NAME}'")) %>%
    group_by(MOTU) %>%
    summarise(across(everything(), ~sum(.x, na.rm = TRUE)))
esv_table_sum <- esv_table_sum %>%
    mutate(total_reads = rowSums(across(starts_with("'${EXPERIMENT}'_'${SAMPLE_NAME}'")))) %>%
    arrange(desc(total_reads)) %>%
    select(-total_reads)
write.table(esv_table_sum, "'${EXPERIMENT}'_'${SAMPLE_NAME}'_motu_table.csv", sep = "\t", row.names = FALSE, quote = FALSE)
'


################################################################################
################################################################################
# Step 6: Taxonomic Assignment

LIB_PREFIX="ULO1"
EXPERIMENT="ULOY"
CORES=5
ORIGINAL_SAMPLE_NAME="10_25_M1_A_C"
SAMPLE_NAME="ULO1_sample_006"

perl ../SOFT/mkLTG/scripts/mkLTG.pl \
    -in ${EXPERIMENT}_${SAMPLE_NAME}_esv_table.csv \
    -taxonomy ../vtam_WOIA_Acari/COInr_WOIAA_taxonomy.tsv \
    -blast_db ../vtam_WOIA_Acari/COInr_WOIAA \
    -outdir out_taxo -out_name ${EXPERIMENT}_${SAMPLE_NAME}_esv_taxo \
    -ltg_params ../vtam_WOIA_Acari/params_truqui.tsv

mv out_taxo/${EXPERIMENT}_${SAMPLE_NAME}_esv_taxo_ltg.tsv .

################################################################################
################################################################################
# Step 7: post-clustering filtering with LULU/mumu
LIB_PREFIX="ULO1"
EXPERIMENT="ULOY"
CORES=5
ORIGINAL_SAMPLE_NAME="10_25_M1_A_C"
SAMPLE_NAME="ULO1_sample_006"

vsearch --usearch_global "${EXPERIMENT}_${SAMPLE_NAME}_otus.fasta"\
    --db "${EXPERIMENT}_${SAMPLE_NAME}_otus.fasta" \
    --self \
    --id 0.84 \
    --iddef 1 \
    --userout "${EXPERIMENT}_${SAMPLE_NAME}_match_list.txt" \
    -userfields query+target+id \
    --maxaccepts 0 \
    --query_cov .9 \
    --maxhits 10

mumu --otu_table "${EXPERIMENT}_${SAMPLE_NAME}_motu_table.csv" \
    --match_list "${EXPERIMENT}_${SAMPLE_NAME}_match_list.txt" \
    --log "${EXPERIMENT}_${SAMPLE_NAME}_mumu.log" \
    --new_otu_table "${EXPERIMENT}_${SAMPLE_NAME}_mumu.tsv"