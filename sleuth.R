library("cowplot")
library("sleuth")
library("dplyr")
library("ggfortify")
library("data.table")

options(scipen = 999) # disable scientific notation

base = "/Users/restifo/NIH/data/kallisto_quant"
metadata = read.table(paste0(base, "/metadata.tsv"), sep = "\t", header = TRUE, stringsAsFactors = FALSE)

# mart = biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl")
# ttg = biomaRt::getBM(attributes = c("ensembl_transcript_id", "transcript_version",
#                                      "ensembl_gene_id", "external_gene_name", "description",
#                                      "transcript_biotype"), mart = mart)
# ttg = dplyr::rename(ttg, target_id = ensembl_transcript_id,
#                      ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

# Load saved biomart index
ttg = readRDS("ttg.Rda")
ttg = dplyr::select(ttg, c("target_id", "ens_gene", "ext_gene"))

# TODO: This sample is corrupted.
metadata = subset(metadata, sample != "4295-BMet-Frag21_FrTu_February_23_2018")

if (!exists("so")) {
  message("creating sleuth object")
  so = sleuth_prep(metadata, target_mapping = ttg, aggregation_column = "ens_gene", extra_bootstrap_summary = TRUE)
  # so = readRDS("sleuth.Rda")
}

# Pre-processing steps
message("getting counts matrix")
matr = sleuth_to_matrix(so, "obs_norm", "est_counts")
matr = t(matr) 
matr = matr[, colSums(matr) != 0] # filter out all columns with column sums of 0. These columns will screw up sleuth's scaling

# Set up PCA
# TODO: can't use autoplot() in code. Must run in the console to actually see the plot.
# ^ this is fixed with wrapping the ggplot() call with plot(). See: https://yihui.name/en/2017/06/top-level-r-expressions/
# m = counts matrix, meta = metadata, colorby = variable to color by
pca = function(m, meta, colorby) {
  if (missing(m) || missing(meta) || missing(colorby)) {
    stop("pca(): need to provide counts matrix, metadata, and a var to color by")
  }
  
  this_pc = prcomp(m, scale. = TRUE, center = TRUE)
  
  the_dt = this_pc$x
  the_names = rownames(the_dt)
  the_dt = as.data.table(the_dt)
  the_dt[, Sample := the_names]
  the_dt = merge(the_dt, meta, by.x = "Sample", by.y = "sample")
  
  plot(ggplot(data = the_dt, aes_string(x = "PC1", y = "PC2", color = colorby)) + geom_point())
}

# "naive" analysis.


