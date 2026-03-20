# This script runs the mjolnir pipeline for the TENE experiment
# this is a merge of the Three different experiments TEN1, TEN2 and TEN3
# I copy the denoised and filtered fasta files and then I keep running from
# ODIN step witht the run_dnoise = FALSE option

# When copying the files from the TEN1, TEN2 and TEN3 experiments,
# make sure to rename them from:
# FTE1_sample_002_ODIN_ESV_Fil1.fasta  to
# TENE_FTE1_sample_002_ODIN_ESV_Fil1.fasta

setwd("data")

library(mjolnir)
experiment <- "ULOY"
cores <- 5
lib_prefix_1 <- c("ULO1", "ULO2", "ULO3", "ULO4")
R1_filenames_1 <- c("ULO1_R1.fastq.gz", "ULO2_R1.fastq.gz", "ULO3_R1.fastq.gz", "ULO4_R1.fastq.gz")
R1_motif = "_R1"
R2_motif = "_R2"

mjolnir1_RAN(experiment = experiment,
            cores = cores,
            lib_prefix = lib_prefix_1,
            R1_filenames = R1_filenames_1,
            R1_motif = R1_motif,
            R2_motif = R2_motif,
            keep_intermediates = TRUE)


mjolnir2_FREYJA(experiment = experiment,
                cores = cores,
                R1_motif = R1_motif,
                R2_motif = R2_motif)


mjolnir3_HELA(experiment = experiment,
               cores = cores)


# experiment = experiment
# cores = cores 
# d = 13
# min_reads_MOTU = 2
# min_reads_ESV = 2
# min_relative = 1 / 50000
# blank_relative = 0.1
# blank_col = "BLANK"
# blank_tag = "blank"
# alpha = 4
# entropy = c(0.47, 0.23, 1.02, 313)
# algorithm = "dnoise_swarm"
# run_dnoise = TRUE
# metadata_table = ''

mjolnir4_ODIN(experiment = experiment, cores = cores, d = 13,
                          min_reads_MOTU = 2, min_reads_ESV = 2,
                          min_relative = 1 / 50000,
                          blank_relative = 0.1,
                          blank_col = "BLANK", blank_tag = "blank",
                          alpha = 4,
                          entropy = c(0.47, 0.23, 1.02, 313),
                          algorithm = "DnoisE_SWARM")



# perl ~/mkLTG/scripts/mkLTG.pl -in TENE_ODIN_counts.tsv -taxonomy ~/mkLTG/vtam_WOIA_Acari/COInr_WOIAA_taxonomy.tsv -blast_db ~/mkLTG/vtam_WOIA_Acari/COInr_WOIAA -outdir out_taxo -out_name TENE_ODIN_counts_taxo -ltg_params ~/mkLTG/params/params_Xavier_truqui.tsv -num_threads 5

# mv out_taxo/TENE_ODIN_counts_taxo_ltg.tsv .

system(paste0("perl ../SOFT/mkLTG/scripts/mkLTG.pl ",
              " -in ", experiment, "_ODIN_counts.tsv ",
              " -taxonomy ../vtam_WOIA_Acari/COInr_WOIAA_taxonomy.tsv ",
              " -blast_db ../vtam_WOIA_Acari/COInr_WOIAA ",
              " -outdir out_taxo ",
              " -out_name ", experiment, "_ODIN_counts_taxo ",
              " -ltg_params ../vtam_WOIA_Acari/params_truqui.tsv ",
              " -num_threads ", cores,
              " ; mv out_taxo/", experiment, "_ODIN_counts_taxo_ltg.tsv ", 
              experiment, "_THOR_annotated.tsv ",
              " ; sed -e 's/^ID/id/g' ", experiment, "_THOR_annotated.tsv > ",
              experiment, "_FRIGGA.tsv ; ",
              " sed -i  's/\\t\\{12\\}/&\\t/g' ", experiment, "_FRIGGA.tsv ; ",
              " sed -i 's/^ID/id/g' ", experiment, "_ODIN_ESV.tsv "
              ), wait = TRUE, intern = FALSE)


mjolnir7_LOKI(experiment = experiment, min_id = .84)

# esv_table <- read.table('TENE_ODIN_ESV.tsv', sep = '\t', header = T) 
# names(esv_table)[names(esv_table) == 'ID'] <- 'id'
# write.table(esv_table, 'TENE_ODIN_ESV.tsv', sep = '\t', quote = F, row.names = F)


metadata_table = ""
output_file = ""
output_file_ESV = ""
min_reads = 0
remove_bacteria = T
remove_contamination = F
contamination_file = "contaminants.txt"
ESV_within_MOTU = T
remove_numts = T

mjolnir8_RAGNAROC(experiment = experiment, remove_numts = T, cores = cores, 
                  remove_contamination = F)


################
# Summary
################

after_ODIN <- read.table(paste0(experiment, "_ODIN_ESV.tsv"),
                        sep = "\t", header = TRUE)

# number of ESVs
num_esvs <- nrow(after_ODIN)
# > num_esvs
# [1] 273635
#number of MOTUs
num_motus <- length(unique(after_ODIN$MOTU))
# > num_motus
# [1] 41076
reads_after_ODIN <- data.frame(reads = colSums(after_ODIN[,grep("sample", names(after_ODIN))]))
