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
