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
  make_option(c("--outdir"), action = 'store', type = 'character', help = '[Required] Output directory for results')
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

input_dir <- paste0("/nfs/turbo/umms-scjp-pank/5_DEG/results/ruvseq_analysis/20250516_ncells", ncells, "_minreads", minreads, "_minprop", minprop, "_withoutAAB/") # change this dir to path to where data is stored on your system

# explore base DESeq2 model
celltype_de_explore <- readRDS(paste0(input_dir, cell.type, "_celltype_de_explore.Rds"))
coldata <- data.frame(colData(celltype_de_explore$`~ diabetes_status_description+sex+age+bmi+ethnicity`$dds))

tmp_df <- distinct(cell.prop[, c("Var1", "Freq")])
colnames(tmp_df) <- c("sample_id", "cell_counts")
coldata <- inner_join(coldata, tmp_df, by = c("sample_id" = "sample_id"))
rownames(coldata) <- coldata$sample_id


# RUVSeq
celltype_ruvseq <- readRDS(paste0(input_dir, cell.type, "_celltype_ruvseq.Rds"))

### corr of just W's on the x-axis
k <- length(celltype_ruvseq)
tmp <- data.frame(celltype_ruvseq[[k]]$W) %>% tibble::rownames_to_column("sample_id")
tmp <- combine_by_sampleid(coldata, tmp)
tmp <- dplyr::select(tmp, -any_of(c("sample_id", "rrid", "samples"))) %>%
    DataExplorer::dummify() %>%
    dplyr::select_if(is_almost_ok, ~.x)
a <- psych::corr.test(tmp[, grep("W_", colnames(tmp))], tmp[, grep("W_|PC", colnames(tmp), invert = T)],
                     use = "na.or.complete", method = "spearman", adjust = "BH")
png(paste0(outdir, cell.type, "_spearmanCorr_onlyWs_vars.png"), width = 15, height = 15, res = 300, units = 'in')
corrplot::corrplot(a$r, p.mat = a$p.adj, tl.cex = .75, sig.level = c(0.01, 0.05, 0.1), pch.cex = 0.8,
                     insig = "label_sig", diag = FALSE)
dev.off()

## calculate which k to search
tmp <- rowSums(a$p.adj < 0.05)
tmp <- data.frame(tmp)
tmp$k <- 1:nrow(tmp)
colnames(tmp) <- c("n_known_vars", "k")
k_stop <- tmp[tmp$n_known_vars == 0, "k"][1] - 1 # get the last k that has significant corr. with known vars
png(paste0(outdir, cell.type, "_nKnownVarsSignificant.png"), width = 8, height = 4, res = 300, units = 'in')
p <- ggplot(tmp, aes(k, n_known_vars)) + geom_line() +
  geom_point() + ggrepel::geom_text_repel(aes(label = glue("{n_known_vars}")), max.overlaps = Inf, box.padding = 0.1) +
  scale_x_continuous(breaks=seq(0, length(celltype_ruvseq), by=1)) +
  labs(x = 'Num. of latent variables', y = 'Total n known vars significantly corr.',
       title = glue("{cell.type}_nCells-{20}_reads-{10}_in-{1/4}-nSamples_fdr-{0.05}"))
print(p)
dev.off()

#find which known vars. corr. with latent vars
tmp <- data.frame(a$p.adj < 0.05)[1:k_stop,]
x <- vector(mode = "list", length = length(1:k_stop))
for (k in 1:k_stop) {
    j <- which(1:k_stop == k)
    x[[j]] <- colnames(tmp)[which(tmp[k,] == TRUE)]
}
names(x) <- paste0("k_", 1:k_stop)
library(UpSetR)
png(paste0(outdir, cell.type, "_overlapKnownVars.png"), width = 8, height = 6, res = 300, units = 'in')
p <- upset(fromList(x), order.by = "freq", nsets = 40,
      text.scale = c(1, 1, 1, 1, 1, 2)) #c(intersection size title, intersection size tick labels,
                                        #set size title, set size tick labels, set names, numbers above bars).
print(p)
dev.off()

#choose k
# Function to find the smallest k with a unique element
find_smallest_k_with_unique <- function(lists) {
  k_with_unique <- c()
  all_elements <- unlist(lists)
  unique_elements <- unique(all_elements)
  element_count <- table(all_elements)
  
  for (i in seq_along(lists)) {
    k_elements <- lists[[i]]
    # Find unique elements in the current list
    unique_in_k <- k_elements[k_elements %in% unique_elements[element_count[k_elements] == 1]]
    if (length(unique_in_k) > 0) {
      k_with_unique <- c(k_with_unique, i)
    }
  }
    if (length(k_with_unique) > 0) {
        return(max(k_with_unique))
    } else {
      return(NULL) # Return NULL if no unique element is found
        }
}

if ( length(grep("diabetes", x)) > 0 ) {
    k_chosen <- grep("diabetes", x) - 1
} else {
    k_chosen <- find_smallest_k_with_unique(x)
}
print(k_chosen)



#### go and kegg functions
suppressPackageStartupMessages(library("clusterProfiler"))
library("org.Hs.eg.db")

### GSEA
for (k in k_chosen) {
    tmp <- celltype_ruvseq[[k]]$de$res[!is.na(celltype_ruvseq[[k]]$de$res$padj),]
    tmp <- tmp[order(tmp$log2FoldChange, decreasing = T),]
    genes <- tmp$log2FoldChange
    names(genes) <- rownames(tmp)
    kruv <- gseGo(genes)
    saveRDS(kruv, paste0(outdir, cell.type, "_k", k, "_gseGO.Rds"))
    kruv <- gseKegg(genes)
    saveRDS(kruv, paste0(outdir, cell.type, "_k", k, "_gseKEGG.Rds"))
}

# choosen k RUV
########## volcano plot
for (k in k_chosen) {
	de_object <- celltype_ruvseq[[k]]$de
	fdr <- 0.05
	xs <- data.frame(de_object$res)
	topxs <- tibble::rownames_to_column(xs[which(xs$padj < fdr), ], var = "geneid")
	t <- unlist(strsplit(x = de_object$design, "+", fixed = T))

	plot <- ggplot(xs, aes(log2FoldChange, -log10(pvalue))) +
    		geom_point(aes(col = ifelse(padj < fdr, "Signif.", "N.S")), size = .5) +
    		scale_color_manual(values = c("gray", "firebrick")) +
    		labs(col = "", title = t[length(t)]) +
    		theme(plot.title = element_text(size = 8))
	plot <- plot + ggrepel::geom_text_repel(data = topxs, aes(x = log2FoldChange, y = -log10(pvalue), label = geneid), size = 3)

	png(paste0(outdir, cell.type, "_k-", k, "_volcanoPlot.png"), width = 6, height = 5, res = 300, units = 'in')
	print(plot)
	dev.off()
}

k <- k_chosen
write.table(celltype_ruvseq[[k]]$de$res, paste0(outdir, cell.type, "_k-chosen_fdr0.05_deg.txt"), row.names = F, sep = "\t", quote = F)
print(cell.type)
print(k_chosen)
d <- many_de_summary(lapply(celltype_ruvseq, function(x) x$de), fdr = 0.05)
print(d[k_chosen,])

