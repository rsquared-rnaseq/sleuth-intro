library(cowplot)
library(sleuth)
library(dplyr)
library(ggfortify)
library(data.table)

options(scipen = 999) # disable scientific notation

metadata_f <- "/data/Robinson-SB/metadata-trials16and14.tsv"
metadata <- read.table(metadata_f, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# mart <- biomaRt::useMart(biomart <- "ENSEMBL_MART_ENSEMBL", dataset <- "hsapiens_gene_ensembl")
# ttg <- biomaRt::getBM(attributes <- c("ensembl_transcript_id", "transcript_version",
#                                      "ensembl_gene_id", "external_gene_name", "description",
#                                      "transcript_biotype"), mart <- mart)
# ttg <- dplyr::rename(ttg, target_id <- ensembl_transcript_id,
#                      ens_gene <- ensembl_gene_id, ext_gene <- external_gene_name)

# Load saved biomart index
ttg <- readRDS("ttg.Rda")
ttg <- dplyr::select(ttg, c("target_id", "ens_gene", "ext_gene"))

# # TODO: These samples are corrupted.
# metadata <- subset(metadata, !(sample %in% c("4295-BMet-Frag21_FrTu_February_23_2018", "4224-1AMet_FrTu_January_20_2017", "4224-2Met_FrTu_January_20_2017", "4267-AMet_FrTu_August_29_2018",
#                                         "4267-BMet_FrTu_August_29_2018", "4267-AMet_FrTu_July_11_2017", "4267-BMet_FrTu_July_11_2017")))

if (!exists("so")) {
  message("creating sleuth object")
  so <- sleuth_prep(metadata, target_mapping = ttg, aggregation_column <- "ens_gene", extra_bootstrap_summary = TRUE)
  # so <- readRDS("sleuth.Rda")
}

# Pre-processing steps
message("getting counts matrix")
matr <- sleuth_to_matrix(so, "obs_norm", "est_counts")
matr <- t(matr) 
matr <- matr[, colSums(matr) != 0] # filter out all columns with column sums of 0. These columns will screw up sleuth's scaling

# Set up PCA
# TODO: can't use autoplot() in code. Must run in the console to actually see the plot.
# ^ this is fixed with wrapping the ggplot() call with plot(). See: https://yihui.name/en/2017/06/top-level-r-expressions/
# m <- counts matrix, meta <- metadata, colorby <- variable to color by
pca <- function(m, meta, colorby) {
  if (missing(m) || missing(meta) || missing(colorby)) {
    stop("pca(): need to provide counts matrix, metadata, and a var to color by")
  }
  
  this_pc <- prcomp(m, scale. = TRUE, center = TRUE)
  
  the_dt <- this_pc$x
  the_names <- rownames(the_dt)
  the_dt <- as.data.table(the_dt)
  the_dt[, Sample := the_names]
  the_dt <- merge(the_dt, meta, by.x = "Sample", by.y = "sample")
  
  ggplot(data = the_dt, aes_string(x = "PC1", y = "PC2", color = colorby)) + geom_point()
}

# "naive" analysis.
# TODO: Identify other covariates to control for
so <- sleuth_fit(so, ~method + sex + age, "reduced")
so <- sleuth_fit(so, ~method + sex + age + response, "full")
# 
# so <- sleuth_lrt(so, "reduced", "full")
# 
# de_genes <- sleuth_results(so, "reduced:full", "lrt", show_all = FALSE)
# de_genes_f <- dplyr::filter(de_genes, qval <= 0.05)
# 
# de_genes_tx <- sleuth_results(so, "reduced:full", "lrt", show_all = FALSE, pval_aggregate = FALSE)
# de_genes_tx_f <- dplyr::filter(de_genes_tx, qval <= 0.05)



