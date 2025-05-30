#' Helper functions, adapted from Walker, Saunders & Rai et al. https://www.nature.com/articles/s41586-023-06693-2

suppressPackageStartupMessages({
  require(magrittr)
  require(ggplot2)
  require(org.Hs.eg.db)
  require(glue)
})

is_almost_ok <- function(vec, tol = .5) {
  (sum(is.na(vec)) <= tol * length(vec)) & (length(unique(vec)) > 1) & (sd(vec) > 0) & (sum(vec != 0) > 1)
}

filter_eval <- function(df, filter_string) {
  dplyr::filter(df, !!friendlyeval::treat_string_as_expr(filter_string))
}

wrap_plot_multi <- function(objs, FUN, ...) {
  lapply(objs, function(x) FUN(x, ...))
}

ep_str_wrap <- function(string, width) {
  x <- gregexpr(' ', string)
  vapply(seq_along(x),
         FUN = function(i) {
           y <- x[[i]]
           n <- nchar(string[i])
           len <- (c(y,n) - c(0, y)) ## length + 1
           idx <- len > width
           j <- which(!idx)
           if (length(j) && max(j) == length(len)) {
             j <- j[-length(j)]
           }
           if (length(j)) {
             idx[j] <- len[j] + len[j+1] > width
           }
           idx <- idx[-length(idx)] ## length - 1
           start <- c(1, y[idx] + 1)
           end <- c(y[idx] - 1, n)
           words <- substring(string[i], start, end)
           paste0(words, collapse="\n")
         },
         FUN.VALUE = character(1)
  )
}

many_de_heatmap <- function(de_objects, fdr = 0.05) {
  print(glue("Using FDR threshold of: {fdr}"))
  tmp <- do.call(rbind, lapply(de_objects, function(x) {
    tmp <- x$res[which(x$res$padj < fdr), ]
    data.frame(design = x$design, gene = rownames(tmp), pvalue = tmp$pvalue, log2FC = tmp$log2FoldChange)
  }))
  rownames(tmp) <- c()
  return(tmp)
}

get_mapping <- function(keys, keytype = "ENSEMBL", column = "SYMBOL") {
  suppressMessages({
    df <- stack(AnnotationDbi::mapIds(
      org.Hs.eg.db, keys = keys, keytype = keytype, column = column, multiVals = "first")
    )
  })
  colnames(df) <- c(column, keytype)
  return(df)
}

#' DESeq2's Variance Stabilizing Transformation
normalize_deseq <- function(counts, ...) {
  # Outputs a "log-like" transformed value
  DESeq2::varianceStabilizingTransformation(counts, ...)
}



get_correlation <- function(df, method = "spearman") {
  dplyr::select(df, -any_of(c("sample_id", "rrid", "samples"))) %>%
    DataExplorer::dummify() %>%
    dplyr::select_if(is_almost_ok, ~.x) %>%
    psych::corr.test(., use = "na.or.complete", method = method, adjust="BH")
}

plot_correlation <- function(df, label, absolute = FALSE) {
  if (absolute) {
    df <- abs(df)
  }
  ComplexHeatmap::Heatmap(df, name = label, border = T, show_column_dend = F, cluster_rows = T, cluster_columns = T,
    row_names_gp = grid::gpar(fontsize = 8), column_names_gp = grid::gpar(fontsize = 8))
}

combine_by_sampleid <- function(df, ...) {
  dplyr::left_join(df, ..., by = c("sample_id"))
}

#' Create DDS object.
#' This object rarely changes, so do not compute it everytime.
create_dds_obj <- function(counts_mat, metadata) {
  sample_ids <- colnames(counts_mat)
  message(glue("info: counts matrix: {paste(dim(counts_mat), collapse = ', ')}"))
  message(glue("info: covariate matrix: {paste(dim(metadata), collapse = ', ')}"))

  col_data <- metadata


  if (any(rownames(col_data) != colnames(counts_mat))) {
    print("names don't match. check or error!")
    return()
  }

  message("info: creating dds object..")
  DESeq2::DESeqDataSetFromMatrix(countData = counts_mat, colData = col_data, design = ~1)
}


run_many_designs_deseq <- function(dds, design, additional_covs, contrast, ...) {
  result <- list()
  for(i in 0:length(additional_covs)) {
    design_str <- paste(c(design, additional_covs[i]), collapse = "+")
    suppressMessages({
      resDe <- run_differential(
        dds,
        design = design_str,
        contrast = contrast_vec,
        ...
      )
    })
    result[[design_str]] <- resDe
  }
  return(result)
}

#' Run the differential expression logic for the given design formula
run_differential <- function(dds, design, contrast, interaction = NULL, fdr = 0.05, shrink = FALSE) {
  if (!is.null(interaction))
    design <- paste0(design, "+", interaction)

  writeLines(glue("info: running design {design}"))
  design(dds) <- as.formula(design)  # replace / update design

  dds <- DESeq2::DESeq(dds)
  vst <- SummarizedExperiment::assay(DESeq2::varianceStabilizingTransformation(dds, blind = FALSE))

  # extract results for the specified FDR threshold
  result <- DESeq2::results(dds, contrast = contrast, alpha = fdr)

  # if shrink is specifies, then report shrunken log2FC using apeglm method
  if (shrink) {
    writeLines("info: estimating shrunk log2FC")
    coef <- gsub(' ', '.', glue("{contrast[1]}_{contrast[2]}_vs_{contrast[3]}"))
    result <- DESeq2::lfcShrink(dds, coef = coef, type = "apeglm")
  }

  return(list(dds = dds, vst = vst, res = result, design = design, shrink = shrink))
}


run_ruvseq <- function(dds, design, contrast, k = 4, p.val.thresh = 0.5, method = "ruvg", ...) {
  # get counts from DESeq2 -- don't normalize yet
  counts <- DESeq2::counts(dds, normalized = FALSE)

  set <- EDASeq::newSeqExpressionSet(counts)
  set <- EDASeq::betweenLaneNormalization(set, which="upper")

  if (method == "ruvg") {
    # Run first-pass DESeq scan to get empirical control genes
    suppressMessages({
      de_res <- run_differential(dds, design = design, contrast = contrast, ...)
    })

    # Use non-significant genes as as empirical control genes
    not.sig <- rownames(de_res$res)[which(de_res$res$pvalue > p.val.thresh)]
    empirical <- rownames(set)[rownames(set) %in% not.sig]
  } else if (method == "ruvs") {
    empirical <- rownames(set)
  }

  print(glue("info: using {length(empirical)} control genes to estimate variation.."))

  result <- list()
  for (num_latent in 1:k) {
    # Run specified method of unwanted variation control
    if (method == "ruvg") {
      print(glue("info: estimating k={num_latent} factors using RUVg.."))
      new_set <- RUVSeq::RUVg(set, empirical, k = num_latent) # Run RUVSeq
    } else if (method == "ruvs") {
      print(glue("info: estimating k={num_latent} factors using RUVs.."))
      groups <- RUVSeq::makeGroups(SummarizedExperiment::colData(dds)[, contrast[1]])
      new_set <- RUVSeq::RUVs(set, empirical, k = num_latent, groups) # Run RUVSeq
    }

    ddsruv <- dds
    for (w in 1:num_latent) {
      col_idx <- glue("W_{w}")
      ddsruv[[col_idx]] <- new_set[[col_idx]]
    }

    new_design <- glue("{design} + {paste('W', 1:num_latent, collapse = '+', sep = '_')}")

    suppressMessages({
      de_res <- run_differential(ddsruv, design = new_design, contrast = contrast, ...)
    })

    result[[new_design]] <- list(de = de_res, W = pData(new_set), normCounts = new_set@assayData$normalizedCounts)
  }

  return(result)
}

plot_ruv_diagnostics <- function(normCounts, metadata, group = "diabetes_status_description", main = NULL, labels = FALSE, ...) {
  ord_idx <- match(colnames(normCounts), metadata$sample_id)
  colors <- RColorBrewer::brewer.pal(9, "Set1")[as.factor(metadata[ord_idx, ][[group]])]
  par(mfrow=c(1,3))
  EDASeq::plotRLE(normCounts, col = colors, outline = FALSE, las = 3, cex.axis = .5, ylab = "Relative Log Expression", main = main, cex.main = 1)
  EDASeq::plotPCA(normCounts, col = colors, cex = 1, cex.lab = 1, cex.axis = 1, labels = labels, main = group, ...)
  
  group = "source"
  colors <- RColorBrewer::brewer.pal(9, "Set1")[as.factor(metadata[ord_idx, ][[group]])]
  EDASeq::plotPCA(normCounts, col = colors, cex = 2, cex.lab = 1, cex.axis = 1, labels = F, main = group)
}

plot_deseq_heatmap <- function(de_object, metadata, split_by = NA, fdr = 0.05, ...) {
  sig.genes <- which(de_object$res$padj < fdr)
  norm_counts <- as.data.frame(de_object$vst[sig.genes, ])

  norm_counts$Gene <- rownames(norm_counts)

  if (nrow(norm_counts) < 1) {
    print("No genes to plot. Adjust FDR if neccessary.")
    return()
  }

  meta_cols <- na.omit(c(split_by, "sample_id", "sex", "chemistry", "study"))

  norm_counts %>%
    dplyr::filter(!is.na(Gene)) %>%
    tidyr::pivot_longer(cols = -c(Gene), names_to = "sample_id", values_to = "Expression") %>%
    dplyr::left_join(., dplyr::select(metadata, any_of(meta_cols)), by = "sample_id") %>%
    { if (!is.na(split_by)) dplyr::group_by(., !!!friendlyeval::treat_strings_as_exprs(split_by)) else . } %>%
    tidyHeatmap::heatmap(
      Gene,
      sample_id,
      Expression,
      .scale = "row",
      palette_value = c("blue", "white", "red"),
      ...
    ) %>%
    tidyHeatmap::add_tile(sex) %>%
    tidyHeatmap::add_tile(chemistry) %>%
    tidyHeatmap::add_tile(study)
}


many_de_summary <- function(de_objects, fdr = 0.05, log2fc.abs = 0) {
  print(glue("Using FDR={fdr} and log2FC={log2fc.abs} to filter genes"))
  do.call(rbind, lapply(de_objects, function(x) {
    tmp <- x$res[which(x$res$padj < fdr & abs(x$res$log2FoldChange) > log2fc.abs), ]
    if (purrr::is_empty(tmp)) {
      return(
        data.frame(total_signif = 0, up = 0, down = 0)
      )
    }
    data.frame(total_signif = nrow(tmp), up = sum(tmp$log2FoldChange > 0), down = sum(tmp$log2FoldChange < 0))
  }))
}


plot_ruv_loadings <- function(obj) {
  tmp_df <- obj$W
  rownames(tmp_df) <- colnames(obj$de$vst)
  ComplexHeatmap::Heatmap(tmp_df, name = " ", cluster_columns = F,
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid::grid.text(sprintf("%.3f", tmp_df[i, j]), x, y, gp = grid::gpar(fontsize = 10))
  })
}

#### go and kegg functions
suppressPackageStartupMessages(library("clusterProfiler"))
library("org.Hs.eg.db")

go <- function(genes) {
  GO <- enrichGO(gene = genes, OrgDb=org.Hs.eg.db, ont = "BP",
                 pvalueCutoff = 1, qvalueCutoff = 1, readable = F, keyType = "SYMBOL")
  res <- as.data.frame(GO@result)
  res$Rank <- seq(1, nrow(res))
  res$GeneRatio <- sapply(strsplit(res$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  res$BgRatio <- sapply(strsplit(res$BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
  res$Fold <- as.numeric(res$GeneRatio)/as.numeric(res$BgRatio)

  return(res)
}

kegg <- function(genes) {
    hs <- org.Hs.eg.db
    my.symbols <- genes
    g <- AnnotationDbi::select(hs, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
    res_kegg <- enrichKEGG(gene = g$ENTREZID, organism = 'hsa',
                           pvalueCutoff = 1, qvalueCutoff = 1)
    res <- as.data.frame(res_kegg@result)
    res$Rank <- seq(1, nrow(res))
    res$GeneRatio <- sapply(strsplit(res$GeneRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
    res$BgRatio <- sapply(strsplit(res$BgRatio, "/"), function(x) as.numeric(x[1])/as.numeric(x[2]))
    res$Fold <- as.numeric(res$GeneRatio)/as.numeric(res$BgRatio)

    return(res)
}

gseKegg <- function(sorted_genes) { # `sorted_genes` has to be sorted
    hs <- org.Hs.eg.db
    my.symbols <- names(sorted_genes)
    g <- AnnotationDbi::select(hs, keys = my.symbols, columns = c("ENTREZID", "SYMBOL"), keytype = "SYMBOL")
    g <- g[!is.na(g$ENTREZID),]
    g <- g[!duplicated(g$SYMBOL),]
    sorted_genes <- sorted_genes[names(sorted_genes) %in% g$SYMBOL]
    names(sorted_genes) <- g$ENTREZID
    res <- gseKEGG(sorted_genes, organism = "hsa", keyType = "kegg", pvalueCutoff = 1)

    return(res)
}

gseGo <- function(sorted_genes) { # `sorted_genes` has to be sorted
    res <- gseGO(sorted_genes, OrgDb=org.Hs.eg.db, ont = "BP", pvalueCutoff = 1, keyType = "SYMBOL")
    return(res)
}

