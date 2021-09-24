# Prepare C084_Hayman_Pasricha data for GEO submission
# Peter Hickey
# 2021-09-24

library(here)
library(SingleCellExperiment)

outdir <- here("GEO")
dir.create(outdir, recursive = TRUE)

# FASTQs -----------------------------------------------------------------------

dir.create(file.path(outdir, "FASTQ"))
# NOTE: These plate-level FASTQ files are created by code/scPipe.R
file.copy(
  from = here("extdata/NN179/scPipe/LCE475/LCE475.R1.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)
file.copy(
  from = here("extdata/NN179/scPipe/LCE475/LCE475.R2.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)
file.copy(
  from = here("extdata/NN183/scPipe/C084_LC475_mini_bulk_back_up_1:10/C084_LC475_mini_bulk_back_up_1:10.R1.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)
file.copy(
  from = here("extdata/NN183/scPipe/C084_LC475_mini_bulk_back_up_1:10/C084_LC475_mini_bulk_back_up_1:10.R2.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)
file.copy(
  from = here("extdata/NN183_NN195/scPipe/LCE_476_BU_8_cycle/LCE_476_BU_8_cycle.R1.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)
file.copy(
  from = here("extdata/NN183_NN195/scPipe/LCE_476_BU_8_cycle/LCE_476_BU_8_cycle.R2.fastq.gz"),
  to = file.path(outdir, "FASTQ"),
  recursive = FALSE,
  overwrite = FALSE)

# SCEs -------------------------------------------------------------------------

dir.create(file.path(outdir, "SCE"))
# NOTE: Starting from the SCE used for the DE analysis (i.e. some samples have
#       been filtered out).
sce <- readRDS(here("data", "SCEs", "C084_Hayman_Pasricha.preprocessed.SCE.rds"))
# NOTE: Restrict data to UMI counts, which were used for the actual analysis.
assays(sce) <- list(counts = assay(sce, "UMI_counts"))
# NOTE: Revert some of the changes to the rowData and colData done in
#       analysis/C084_Hayman_Pasricha.preprocess.Rmd
sce$sex_colours <- NULL
sce$timepoint_colours <- NULL
sce$treatment_colours <- NULL
sce$plate_number_colours <- NULL
sce$batch <- NULL # Redundant with plate_number
rownames(sce) <- rowData(sce)$ENSEMBL.GENEID
rowData(sce) <- S4Vectors::make_zero_col_DFrame(nrow(sce))

# Gene counts
write.csv(
  x = as.data.frame(as.matrix(counts(sce))),
  file = gzfile(file.path(outdir, "SCE", "gene_counts.csv.gz")),
  row.names = TRUE)

# ERCC counts
write.csv(
  x = as.data.frame(as.matrix(counts(altExp(sce)))),
  file = gzfile(file.path(outdir, "SCE", "ERCC_counts.csv.gz")),
  row.names = TRUE)

# colData
write.csv(
  x = as.data.frame(colData(sce)),
  file = gzfile(file.path(outdir, "SCE", "sample_sheet.csv.gz")),
  row.names = TRUE)
