# Helper function to coerce a DataFrame to a data.frame while preserving column
# names as is.
.adf <- function(x) {
  setNames(as.data.frame(x), colnames(x))
}

# Helper function to Combine data from 2 SCEs using gene names.
# NOTE: This assumes more than I'd like about the rowData and doesn't do much
#       checking of these assumptions.
.combine <- function(x, y, rowData_by = c("ENSEMBL", "SYMBOL", "CHR")) {
  if (is.null(rowData_by)) {
    rowData <- dplyr::full_join(
      .adf(rowData(x)) %>%
        tibble::rownames_to_column(var = "gene"),
      .adf(rowData(y)) %>%
        tibble::rownames_to_column(var = "gene")) %>%
      tibble::column_to_rownames("gene") %>%
      DataFrame(., row.names = rownames(.))
  } else {
    rowData <- dplyr::full_join(
      .adf(rowData(x)[, rowData_by, drop = FALSE]),
      .adf(rowData(y)[, rowData_by, drop = FALSE]),
      by = rowData_by) %>%
      DataFrame(row.names = scater::uniquifyFeatureNames(
        .$ENSEMBL,
        .$SYMBOL))
    rownames(x) <- rownames(rowData)[match(rowData(x)$ENSEMBL, rowData$ENSEMBL)]
    rownames(y) <- rownames(rowData)[match(rowData(y)$ENSEMBL, rowData$ENSEMBL)]
  }

  colData <- rbind(colData(x), colData(y))

  counts <- matrix(
    data = 0L,
    nrow = nrow(rowData), ncol = nrow(colData),
    dimnames = list(rownames(rowData), rownames(colData)))
  counts[rownames(x), colnames(x)] <- counts(
    x,
    withDimnames = FALSE)
  counts[rownames(y), colnames(y)] <- counts(
    y,
    withDimnames = FALSE)

  stopifnot(
    identical(
      metadata(x)$scPipe$version,
      metadata(y)$scPipe$version))
  stopifnot(
    identical(
      metadata(x)$scPipe$QC_cols,
      metadata(y)$scPipe$QC_cols))
  stopifnot(
    identical(
      metadata(x)$scPipe$demultiplex_info$status,
      metadata(y)$scPipe$demultiplex_info$status))
  stopifnot(
    identical(
      metadata(x)$scPipe$UMI_dup_info$duplication.number,
      metadata(y)$scPipe$UMI_dup_info$duplication.number))
  stopifnot(identical(metadata(x)$Biomart, metadata(y)$Biomart))
  metadata <- list(
    scPipe = list(
      version = metadata(x)$scPipe$version,
      QC_cols = metadata(x)$scPipe$QC_cols,
      demultiplex_info = data.frame(
        status = metadata(x)$scPipe$demultiplex_info$status,
        count = metadata(x)$scPipe$demultiplex_info$count +
          metadata(y)$scPipe$demultiplex_info$count),
      UMI_dup_info = data.frame(
        duplication.number = metadata(
          x)$scPipe$UMI_dup_info$duplication.number,
        count = metadata(x)$scPipe$UMI_dup_info$count +
          metadata(y)$scPipe$UMI_dup_info$count)),
    Biomart = metadata(x)$Biomart)

  sce <- SingleCellExperiment(
    rowData = rowData,
    colData = colData,
    assays = list(counts = counts),
    metadata = metadata)

  stopifnot(identical(int_metadata(x), int_metadata(y)))
  int_metadata(sce) <- int_metadata(x)
  int_elementMetadata <- dplyr::full_join(
    x = .adf(int_elementMetadata(x)) %>%
      tibble::add_column(gene = rownames(x)),
    y = .adf(int_elementMetadata(y)) %>%
      tibble::add_column(gene = rownames(y))) %>%
    tibble::column_to_rownames("gene") %>%
    DataFrame()
  int_elementMetadata(sce) <- int_elementMetadata

  stopifnot(validObject(sce))
  sce
}

.cbindSCEs <- function(list_of_sce, rowData_by = 1:6) {
  do.call(
    cbind,
    lapply(list_of_sce, function(sce) {
      # NOTE: Some fudging to combine only the necessary bits of each SCE
      #       (basically, don't include any QC metrics).
      rowData(sce) <- rowData(sce)[, rowData_by]
      sce
    }))
}

# NOTE: My best guess of what this function does https://github.com/MarioniLab/compareSingleCell/blob/543aa28e3ae25fad4ffb1d47c27c8a364966095c/vignettes/embryo_expression.Rmd#L76
.sumCountsAcrossCells <- function(sce, cluster_sample) {
  counts <- counts(sce, withDimnames = FALSE)
  edgeR::sumTechReps(counts, cluster_sample)
}

# NOTE: Need to use my own gene counting function because not using UMI
#       deduplication.
geneCountingNoUMIDedup <- function(outdir, bc_anno) {
  files <- list.files(file.path(outdir, "count"), full.names = TRUE)
  names(files) <- sub("\\.csv", "", basename(files))
  counts <- lapply(files, function(file) {
    message(basename(file))
    data.table::fread(file, select = 1)[, table(gene_id)]
  })
  genes <- Reduce(union, lapply(counts, names))
  x <- matrix(
    0L,
    nrow = length(genes),
    ncol = length(files),
    dimnames = list(genes, names(counts)))
  for (j in names(counts)) {
    xx <- counts[[j]]
    x[names(xx), j] <- xx
  }
  z <- cbind(
    data.frame(gene_id = rownames(x)),
    as.data.frame(x))
  data.table::fwrite(
    x = z,
    file = file.path(paste0(outdir, "_no_dedup"), "gene_count.csv"),
    row.names = FALSE,
    nThread = 1)
}

# Take a DataFrame with AtomicList columns and return a DataFrame where these
# columns have been flattened by paste-ing together the elements separated by
# `sep`.
flattenDF <- function(x, sep = "; ") {
  DataFrame(
    endoapply(x, function(xx) {
      if (!is(xx, "AtomicList")) {
        return(xx)
      }
      unstrsplit(as(xx, "CharacterList"), sep = sep)
    }),
    row.names = rownames(x))
}

# NOTE: This function is customised for the C084_Hayman_Pasricha project.
createDEGOutputs <- function(outdir, efit, x, prefix, suffix, o, groups, sample.cols) {

  message("Creating CSVs of DEGs")
  for (label in colnames(efit)) {
    message("\t", label)
    gzout <- gzfile(
      description = file.path(
        outdir,
        paste(c(label, suffix, "DEGs.csv.gz"), collapse = ".")),
      open = "wb")
    write.csv(
      topTable(efit, coef = label, n = Inf),
      gzout,
      # NOTE: quote = TRUE needed because some fields contain commas.
      quote = TRUE,
      row.names = FALSE)
    close(gzout)
  }

  # NOTE: Doesn't make sense to create Venn diagram if there's only one
  #       contrasts
  if (ncol(efit) > 1) {
    message("Creating Venn diagram")
    dt <- decideTests(efit)
    par(mfrow = c(1, 1))
    pdf(
      here(
        "output",
        "DEGs",
        paste(c(prefix, suffix, "Venn_diagram.pdf"), collapse = ".")),
      width = 7,
      height = 7)
    vennDiagram(dt, cex = 0.8)
    dev.off()
  } else {
    message("Skipping Venn diagram")
  }

  message("Creating heatmaps")
  for (label in colnames(efit)) {
    all_features <- rownames(
      topTable(efit, coef = label, p.value = 0.05, n = Inf))
    # Filter genes
    features <- setdiff(
      all_features,
      c(mito_set, ribo_set, pseudogene_set, sex_set))
    # Select top-100
    features <- head(features, 100)
    if (length(features) > 1) {
      mat <- cpm(x, log = TRUE)[features, ]
      # NOTE: Row-normalize within each plate_number
      mat[, x$samples$plate_number == "LC475"] <-
        mat[, x$samples$plate_number == "LC475"] -
        rowMeans(mat[, x$samples$plate_number == "LC475"])
      mat[, x$samples$plate_number == "LC476"] <-
        mat[, x$samples$plate_number == "LC476"] -
        rowMeans(mat[, x$samples$plate_number == "LC476"])
      pheatmap(
        mat = mat[, o],
        color = hcl.colors(101, "Blue-Red 3"),
        cluster_rows = TRUE,
        cluster_cols = FALSE,
        show_colnames = FALSE,
        annotation_col = data.frame(
          plate_number = x$samples$plate_number[o],
          treatment = x$samples$treatment[o],
          timepoint = x$samples$timepoint[o],
          sex = x$samples$sex[o],
          row.names = colnames(x)[o]),
        annotation_colors = list(
          plate_number = plate_number_colours,
          treatment = treatment_colours,
          timepoint = timepoint_colours,
          sex = sex_colours),
        breaks = seq(-max(abs(mat)), max(abs(mat)), length.out = 101),
        fontsize = 6,
        gaps_col = cumsum(runLength(Rle(x$samples$plate_number[o]))),
        main = label,
        filename = here(
          "output",
          "DEGs",
          paste(c(label, suffix, "heatmap.pdf"), collapse = ".")))
    }
  }

  message("Creating Glimma plots")
  for (label in colnames(efit)) {
    # NOTE: Doesn't make sense to create Glimma plots for numeric covariates.
    if (label %in% c("HAZ24", "WHZ24", "WAZ24", "hb24")) {
      next
    }
    message("\t", label)
    anno <- cbind(
      DataFrame(GeneID = rownames(x)),
      x$genes)
    colnames(anno) <- sub("ENSEMBL", "ENS", colnames(anno))
    colnames(anno) <- sub("\\.GENE", "\\.", colnames(anno))
    Glimma::glMDPlot(
      x = efit[, label],
      counts = x,
      groups = groups,
      status = decideTests(efit)[, label],
      transform = TRUE,
      main = label,
      anno = anno,
      display.columns = c(
        "GeneID",
        paste0("ENS.", c("ID", "BIOTYPE", "SEQNAME")),
        paste0("NCBI.", c("ALIAS", "NAME"))),
      sample.cols = sample.cols,
      path = here("output"),
      html = paste(c(label, suffix, "md-plot"), collapse = "."),
      launch = FALSE)
  }

  message("Creating CSVs of GO and KEGG analyses")
  for (label in colnames(efit)) {
    # NOTE: Features map to multiple ENTREZID, so just select the first.
    geneid <- sapply(strsplit(efit$genes$NCBI.ENTREZID, ";"), "[[", 1)

    # Run GO analysis
    go <- limma::goana(efit, coef = label, species = "Hs", geneid = geneid)

    # Get gene symbols associated with each GO term.
    library(org.Hs.eg.db)
    go_symbol_df <- select(
      org.Hs.eg.db,
      keys = rownames(go),
      columns = "SYMBOL",
      # NOTE: `keytype = "GO"` returns 'incomplete' results whereas
      #       `keytype = "GOALL"` returns the 'complete' results.
      #       I think this is the same issue described in the goseq vignette
      #       (section 6.7):
      #       "It is important to use the org.*.egGO2ALLEGS and NOT the org.*.egGO
      #       object to create the mapping between GO categories and gene
      #        identifiers, as the latter does not include the links to genes
      #        arising from 'child' GO categories."
      keytype = "GOALL")
    go_symbol_df <- dplyr::left_join(
      go_symbol_df,
      as.data.frame(
        topTable(
          efit,
          coef = label,
          n = Inf))[, c("logFC", "adj.P.Val", "NCBI.SYMBOL")],
      by = c("SYMBOL" = "NCBI.SYMBOL"))

    # Extract gene symbols in each GO term that are up- or down-regulated.
    go <- dplyr::left_join(
      go %>% tibble::rownames_to_column("GO"),
      go_symbol_df %>%
        dplyr::group_by(GOALL) %>%
        dplyr::summarise(
          upGeneSymbols =
            paste(
              unique(
                SYMBOL[which(logFC > 0 & adj.P.Val < 0.05)]), collapse = ", "),
          downGeneSymbols =
            paste(
              unique(
                SYMBOL[which(logFC < 0 & adj.P.Val < 0.05)]), collapse = ", ")),
      by = c("GO" = "GOALL"))
    go <- tibble::column_to_rownames(go, "GO")

    # Write output
    gzout <- gzfile(
      description = file.path(
        outdir,
        paste(c(label, suffix, "GO.csv.gz"), collapse = ".")),
      open = "wb")
    write.csv(
      go,
      gzout,
      # NOTE: quote = TRUE needed because some fields contain commas.
      quote = TRUE,
      row.names = TRUE)
    close(gzout)

    # Run KEGG analysis
    kegg <- kegga(efit, coef = label, species = "Hs", geneid = geneid)

    kegg_symbol_df <- select(
      org.Hs.eg.db,
      # NOTE: The org.Hs.eg.db does not use the 'path:hsa'-prefix.
      keys = sub("path:hsa", "", rownames(kegg)),
      columns = "SYMBOL",
      keytype = "PATH")
    kegg_symbol_df <- dplyr::left_join(
      kegg_symbol_df,
      as.data.frame(
        topTable(
          efit,
          coef = label,
          n = Inf))[, c("logFC", "adj.P.Val", "NCBI.SYMBOL")],
      by = c("SYMBOL" = "NCBI.SYMBOL"))

    # Extract gene symbols in each KEGG pathway that are up- or down-regulated.
    kegg <- dplyr::left_join(
      kegg %>%
        tibble::rownames_to_column("PATH") %>%
        dplyr::mutate(PATH = sub("path:hsa", "", PATH)),
      kegg_symbol_df %>%
        dplyr::group_by(PATH) %>%
        dplyr::summarise(
          upGeneSymbols =
            paste(
              unique(
                SYMBOL[which(logFC > 0 & adj.P.Val < 0.05)]), collapse = ", "),
          downGeneSymbols =
            paste(
              unique(
                SYMBOL[which(logFC < 0 & adj.P.Val < 0.05)]), collapse = ", ")),
      by = c("PATH" = "PATH"))
    kegg <- kegg %>%
      dplyr::mutate(PATH = paste0("path:hsa", PATH)) %>%
      tibble::column_to_rownames("PATH")

    # Write output
    gzout <- gzfile(
      description = file.path(
        outdir,
        paste(c(label, suffix, "KEGG.csv.gz"), collapse = ".")),
      open = "wb")
    write.csv(
      kegg,
      gzout,
      # NOTE: quote = TRUE needed because some fields contain commas.
      quote = TRUE,
      row.names = TRUE)
    close(gzout)
  }
}
