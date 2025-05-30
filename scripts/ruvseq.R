suppressMessages(library(Seurat))
suppressMessages(library(stringr))
suppressMessages(library(parallel))
suppressMessages(library(readr))
suppressMessages(library(DESeq2))
suppressMessages(library(limma))
suppressMessages(library(edgeR))
suppressMessages(library(GenomicFeatures))
suppressMessages(library(data.table))
library(ggplot2)
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(S4Vectors))
suppressPackageStartupMessages(library(GenomicRanges))
library(patchwork)
library(biomaRt)

suppressPackageStartupMessages(library("RUVSeq"))
suppressPackageStartupMessages(library("ComplexHeatmap"))
suppressPackageStartupMessages(library("friendlyeval"))

source("/nfs/turbo/umms-scjp-pank/5_DEG/scripts/PanKbase-DEG-analysis/scripts/0_utils.R")
library(optparse)

option_list <- list(
  make_option(c("--celltype"), action = 'store', type = 'character', help = '[Required] Cell type'),
  make_option(c("--ncells"), action = 'store', type = 'numeric', default = 20, help = '[Optional] Min n cells to include a samples, default: 20'),
  make_option(c("--minreads"), action = 'store', type = 'numeric', default = 10, help = '[Optional] Min n reads to keep a gene'),
  make_option(c("--minprop"), action = 'store', type = 'numeric', default = 0.25, help = '[Optional] Min proportion of samples that have minreads reads to keep a gene'),
  make_option(c("--nlatent"), action = 'store', type = 'numeric', default = 30, help = '[Optional] Number of latent vars to test; should be smaller than total number of samples minus number of known vars in model'),
  make_option(c("--outdir"), action = 'store', type = 'character', help = '[Required] Directory to save outputs')
)

option_parser <- OptionParser(usage = "usage: Rscript %prog [options]", option_list = option_list, add_help_option = T, description = '')
opts <- parse_args(option_parser)

cell.type <- opts$celltype
ncells <- opts$ncells
minreads <- opts$minreads
minprop <- opts$minprop
outdir <- opts$outdir
nlatent <- opts$nlatent

samples <- readRDS("/nfs/turbo/umms-scjp-pank/5_DEG/results/samplelist.rds") # sample list
meta_in_sc <- readRDS("/nfs/turbo/umms-scjp-pank/5_DEG/results/metadata_for_DEG.rds")
cell.prop <- read.table(paste0("/nfs/turbo/umms-scjp-pank/5_DEG/results/cell_proportion/cell.prop_", cell.type, ".txt"), header = T) # cell counts and proportions in large map
metadata <- meta_in_sc
metadata <- metadata[metadata$treatments == "no_treatment",]

a <- data.frame(table(metadata$rrid))
coldata <- metadata[metadata$rrid %in% a[a$Freq == 1, "Var1"],]
for (i in a[a$Freq > 1, "Var1"]) {
    tmp <- metadata[metadata$rrid == i,]
    set.seed(1234)
    tmp <- tmp[sample(1:nrow(tmp), 1),]
    coldata <- rbind(coldata, tmp) #coldata has unique sample per donor; no donor has multiple samples
}

cell.prop <- inner_join(cell.prop, metadata, by = c("Var1" = "samples"))
cell.prop$aab <- "0"
cell.prop$aab <- ifelse(cell.prop$aab_gada == "TRUE", "AAB+", cell.prop$aab)
cell.prop$aab <- ifelse(cell.prop$aab_ia_2 == "TRUE", "AAB+", cell.prop$aab)
cell.prop$aab <- ifelse(cell.prop$aab_iaa == "TRUE", "AAB+", cell.prop$aab)
cell.prop$aab <- ifelse(cell.prop$aab_znt8 == "TRUE", "AAB+", cell.prop$aab)
cell.prop$diabetes_status_description <- ifelse(cell.prop$aab == "0" & cell.prop$diabetes_status_description != "T1DM",
                                                cell.prop$diabetes_status_description,
                                                ifelse(cell.prop$aab == "AAB+" & cell.prop$diabetes_status_description == "NonDiabetic",
                                                       cell.prop$aab, cell.prop$diabetes_status_description))
tmp <- cell.prop[cell.prop$Freq > ncells,]
tmp$diabetes_status_description <- factor(tmp$diabetes_status_description, levels = c("NonDiabetic", "PreDiabetic", "AAB+", "T1DM", "T2DM"))
tmp <- cell.prop[cell.prop$Freq > ncells & cell.prop$diabetes_status_description %in% c("NonDiabetic", "T1DM") &
          cell.prop$Var1 %in% intersect(coldata$samples, tmp$Var1),]
coldata <- coldata[coldata$samples %in% tmp$Var1,]
coldata <- data.frame(coldata)
rownames(coldata) <- coldata$samples

# base DESeq2
print("start DESeq2 base exploration")

dir <- '/nfs/turbo/umms-scjp-pank/5_DEG/results/pseudobulk_counts/'
raw_mat <- read.table(paste0(dir, cell.type, "_sample_gex_total_counts.txt"), header = T)
colnames(raw_mat) <- samples
raw_mat <- raw_mat[, intersect(coldata$samples, tmp$Var1)] # keep only samples with > `n_cells` cells
raw_mat <- raw_mat[, which(colSums(raw_mat) > 0)] # remove samples that do not have any cells in the population

# print out the number of samples per diabetes status
## to update to print out more attributes such as age groups, sex
print(table(tmp[tmp$Var1 %in% colnames(raw_mat), "diabetes_status_description"]))

# filter genes that have some min_reads raw counts in at least min_prop ratio of samples
basic_filter <- function (row, min_reads = minreads, min_prop = minprop) {
  mean(row >= min_reads) >= min_prop
}
keep <- apply(raw_mat, 1, basic_filter) #1 means apply the function to each row
raw_mat <- raw_mat[keep, ]
coldata <- coldata[colnames(raw_mat),]
print("dim(coldata)")
print(dim(coldata))
print("dim(raw_mat)")
print(dim(raw_mat))

print("check if all rownames(coldata) is in colnames(raw_mat):")
print(all(rownames(coldata) %in% colnames(raw_mat)))

# check if rownames(coldata) == colnames(raw_mat)
print("check if rownames(coldata) == colnames(raw_mat):")
print(all(rownames(coldata) == colnames(raw_mat)))

for (i in c('study', 'samples', 'rrid', 'treatments', 'chemistry', 'sex', 
            'diabetes_status_description', 'tissue_source', 'ethnicity',
            'aab_gada', 'aab_ia_2', 'aab_iaa', 'aab_znt8')) {
    coldata[, i] <- as.factor(coldata[, i])
    
}
colnames(coldata)[which(colnames(coldata) == "samples")] <- "sample_id"
dds <- DESeqDataSetFromMatrix(countData = raw_mat,
                              colData = coldata)

celltype_raw_counts <- counts(dds)

# normalize using varianceStabilizingTransform from DESeq2, blind = TRUE
celltype_vst_counts <- normalize_deseq(celltype_raw_counts)

# normalize for library sizes
size.factors <- DESeq2::estimateSizeFactorsForMatrix(celltype_raw_counts)
celltype_norm_counts <- t(apply(celltype_raw_counts, 1, function(x) x/size.factors))
celltype_norm_counts.log <- log2(celltype_norm_counts + 1) # add a pseudocount

png(paste0(outdir, cell.type, "_raw__RLE.png"), width = 10, height = 6, res = 300, units = 'in')
ord_idx <- match(colnames(celltype_norm_counts.log), coldata$sample_id)
group = "diabetes_status_description"
par(mfrow=c(2,3))
colors <- RColorBrewer::brewer.pal(9, "Set1")[as.factor(coldata[ord_idx, ][[group]])]
EDASeq::plotRLE(celltype_norm_counts, col = colors, outline = FALSE, las = 3, cex.axis = 1, ylab = "Relative Log Expression", main = NULL, cex.main = .5)
EDASeq::plotPCA(celltype_norm_counts, col = colors, cex = .5, cex.lab = 0.75, cex.axis = 1, labels = F, main = group)
group = "tissue_source"
colors <- RColorBrewer::brewer.pal(9, "Set1")[as.factor(coldata[ord_idx, ][[group]])]
EDASeq::plotPCA(celltype_norm_counts, col = colors, cex = .5, cex.lab = 0.75, cex.axis = 1, labels = T, main = group)
group = "sex"
colors <- RColorBrewer::brewer.pal(9, "Set1")[as.factor(coldata[ord_idx, ][[group]])]
EDASeq::plotPCA(celltype_norm_counts, col = colors, cex = .5, cex.lab = 0.75, cex.axis = 1, labels = F, main = group)
group = "chemistry"
colors <- RColorBrewer::brewer.pal(9, "Set1")[as.factor(coldata[ord_idx, ][[group]])]
EDASeq::plotPCA(celltype_norm_counts, col = colors, cex = .5, cex.lab = 0.75, cex.axis = 1, labels = F, main = group)
dev.off()

dds <- create_dds_obj(celltype_raw_counts, coldata)
base_design <- "~ diabetes_status_description+sex+age+bmi+ethnicity"
additional_covs <- c("chemistry")
contrast_vec = c("diabetes_status_description", "T1DM", "NonDiabetic") # fold change = numerator / denominator
celltype_de_explore <- run_many_designs_deseq(
  dds,
  base_design,
  additional_covs,
  contrast = contrast_vec,
  shrink = FALSE
)
saveRDS(celltype_de_explore, paste0(outdir, cell.type, "_celltype_de_explore.Rds"))

# RUVseq
print("starting RUVSeq")
dds <- create_dds_obj(
  celltype_raw_counts,
  coldata
)
celltype_ruvseq <- run_ruvseq(
    dds,
    design = "~ diabetes_status_description+sex+age+bmi+ethnicity+chemistry",
    contrast = contrast_vec,
    k = min(nlatent, 30), #nrow(coldata) = number of samples - 6 covariates in base model - 3 to be sure
    p.val.thresh = 0.5,
    method = "ruvg"
)
saveRDS(celltype_ruvseq, paste0(outdir, cell.type, "_celltype_ruvseq.Rds"))

