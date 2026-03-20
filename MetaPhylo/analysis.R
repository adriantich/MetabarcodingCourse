# In this script we will run some data analisis with the data from 
# "Metabarcoding reveals high-­resolution biogeographical and
# metaphylogeographical patterns through marine barriers"
# by Antich et al 2022 (DOI: 10.1111/jbi.14548)


# Load the data of the ESV table and the metadata
setwd("MetaPhylo")

library(vegan)
library(dplyr)
library(MetaPhyloTools)
library(ggplot2)

df <- read.csv("ESV.csv")
metadata <- read.csv("metadata.csv")
sample_names <- metadata$codi

df <- df[rowSums(df[, sample_names]) > 0, ] # remove ESVs with zero counts across all samples
df_rarefied <- df
df_rarefied[, sample_names] <- t(rrarefy(t(df_rarefied[, sample_names]), min(colSums(df_rarefied[, sample_names]))))
djost_grouped <- pairwise_djost(
    df = df_rarefied,
    sample_names = sample_names,
    motu_col = "motu",
    rarefy = TRUE,
    rarefy_to = min
)

set.seed(123)
nmds_djost <- metaMDS(djost_grouped[["mean_matrix"]])

nmds_data <- data.frame(
    x = nmds_djost$points[, 1],
    y = nmds_djost$points[, 2],
    sample = rownames(nmds_djost$points),
    codi_comunitat = metadata$codi.comunitat,
    Area_filogenetica = metadata$codi.comunitat
)

color_palette <- hsv(c(0,0,0,0,0.33,0.33,0.33,0.33,0.66,0.66,0.66,0.66), #color
                     c(1,.8,.6,.4,1,.8,.6,.4,1,.8,.6,.4), #saturacio
                     c(.5,1,1,1,.5,1,1,1,.5,1,1,1)) #brillo

order_factor_levels <- function(x){
  if (!is.numeric(x)&!is.logical(x)) {
    x <- factor(x,levels = unique(x))
  } else {
    x
  }
}

nmds_data$sample <- order_factor_levels(nmds_data$sample)
nmds_data$Area_filogenetica <- order_factor_levels(nmds_data$Area_filogenetica)
nmds_data$codi_comunitat <- order_factor_levels(nmds_data$codi_comunitat)


ggplot(nmds_data, aes(x = x, y = y)) +
    geom_point(aes(color = codi_comunitat), size = 8, shape = 16) +
    scale_color_manual(
        "Outside",
        values = color_palette
    ) +
    theme_minimal() +
    theme(
        # legend.box = "horizontal", # Place legends side by side
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
ggsave("nmds_djost.png", width = 8, height = 6)
####################################################################
# join the reads of the rows assigned to the same value in "motu"
df_grouped_motu <- df_rarefied %>%
    group_by(motu) %>%
    summarise(across(all_of(sample_names), sum))

df_grouped_motu[,sample_names] <- t(decostand(t(df_grouped_motu[, sample_names]), method = "hellinger"))

nmds_BC <- metaMDS(t(df_grouped_motu[, sample_names]))

nmds_data_BC <- data.frame(
    x = nmds_BC$points[, 1],
    y = nmds_BC$points[, 2],
    sample = rownames(nmds_BC$points),
    codi_comunitat = metadata$codi.comunitat,
    Area_filogenetica = metadata$codi.comunitat
)

color_palette <- hsv(c(0,0,0,0,0.33,0.33,0.33,0.33,0.66,0.66,0.66,0.66), #color
                     c(1,.8,.6,.4,1,.8,.6,.4,1,.8,.6,.4), #saturacio
                     c(.5,1,1,1,.5,1,1,1,.5,1,1,1)) #brillo

order_factor_levels <- function(x){
  if (!is.numeric(x)&!is.logical(x)) {
    x <- factor(x,levels = unique(x))
  } else {
    x
  }
}
  
nmds_data_BC$sample <- order_factor_levels(nmds_data_BC$sample)
nmds_data_BC$Area_filogenetica <- order_factor_levels(nmds_data_BC$Area_filogenetica)
nmds_data_BC$codi_comunitat <- order_factor_levels(nmds_data_BC$codi_comunitat)


ggplot(nmds_data_BC, aes(x = x, y = y)) +
    geom_point(aes(color = codi_comunitat), size = 8, shape = 16) +
    scale_color_manual(
        "Outside",
        values = color_palette
    ) +
    theme_minimal() +
    theme(
        # legend.box = "horizontal", # Place legends side by side
        axis.text.x = element_text(angle = 45, hjust = 1)
    )
ggsave("nmds_BC.png", width = 8, height = 6)


####################################################################
# Haplotype networks within MOTUs
for (motu in unique(df_rarefied$motu)[1:10]) {
    motu_tab <- df_rarefied[df_rarefied$motu == motu, ]
    if (nrow(motu_tab) == 0) {
        next
    }
    rowSums(motu_tab[,sample_names])

    haploNet_data_object <- haploNet_data(motu_tab,
        metadata = metadata,
        grouping_col = "codi.comunitat",
        sample_col = "codi",
        seq_col = "seq",
        id_col = "id"
    )
    haplot_pegas <- haplodata4ggplot(
        haploNet_data_object = haploNet_data_object,
        method = "pegas"
    )
    haplot_fruch <- haplodata4ggplot(
        haploNet_data_object = haploNet_data_object,
        method = "fruchtermanreingold"
    )
    plot_pegas <- haplo_ggplot(data = haplot_pegas, legend_pos = "none")
    plot_pegas <- plot_pegas + scale_fill_manual(
        values = color_palette
    ) +
        theme_minimal() +
        theme(title = element_text(size = 24))
    ggsave(paste0("networks/haplo_network_pegas_", motu, ".png"), plot = plot_pegas, width = 10, height = 10)
    plot_fruch <- haplo_ggplot(data = haplot_fruch, legend_pos = "none")
    plot_fruch <- plot_fruch + scale_fill_manual(
        values = color_palette
    ) +
        theme_minimal() +
        theme(title = element_text(size = 24))
    ggsave(paste0("networks/haplo_network_fruc_", motu, ".png"), plot = plot_fruch, width = 10, height = 10)
}
