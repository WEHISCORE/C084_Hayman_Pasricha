# Process NN179 and NN183 (C084) with scPipe
# Peter Hickey
# 2020-02-24

# Setup ------------------------------------------------------------------------

library(scPipe)
library(Rsubread)
library(here)
library(readxl)
library(dplyr)
library(janitor)

options("mc.cores" = 1L)

source(here("code", "helper_functions.R"))

# Construct NN179 sample sheet -------------------------------------------------

file_nn179 <- here(
  "data",
  "sample_sheets",
  "C084_Naik_SaraTomei_NN179_SeqPrimer_layout_Jan20.xlsx")

# NOTE: Header row is split across 2 lines, which I combine into 1 before
#       reading in the rest of the spreadsheet.
header_row <- read_excel(
  path = file_nn179,
  sheet = "Sample & Index",
  skip = 2,
  n_max = 1)

header_row <- paste0(colnames(header_row), header_row[1, ])
header_row <- gsub("^\\.\\.\\.[0-9]+", "", header_row)
sample_sheet_nn179 <- read_excel(
  path = file_nn179,
  sheet = "Sample & Index",
  skip = 4,
  col_names = header_row,
  # NOTE: Setting the max guess_max value avoids problems with incorrectly
  #       guessed columns
  #       (https://github.com/tidyverse/readxl/issues/414#issuecomment-352437730)
  guess_max = 1048576)
# Tidy up names
sample_sheet_nn179 <- clean_names(sample_sheet_nn179)
# NOTE: Restrict to the plates for this project.
sample_sheet_nn179 <- sample_sheet_nn179 %>%
  filter(plate_number == "LCE 475")
# Remove empty rows/columns.
sample_sheet_nn179 <- remove_empty(sample_sheet_nn179)

# Filter out those empty wells.
sample_sheet_nn179 <- sample_sheet_nn179 %>%
  filter(sample_name != "empty")
# Drop remaining index sorting column (no index sorting on these samples)
sample_sheet_nn179 <- sample_sheet_nn179 %>%
  select(-indexing_fsc_a)

# Some final tidying.
sample_sheet_nn179 <- sample_sheet_nn179 %>%
  mutate(
    plate_number = sub(" ", "", plate_number),
    # NOTE: There are some wonky well_positions (e.g., 'I19=A1') that need to
    #       be fixed (these occur because it means well I19 with primer A1,
    #       in SCORE's terminology. I've asked for this to be avoided going
    #       forward.).
    well_position = gsub(" ", "", well_position),
    well_position = sapply(strsplit(well_position, "="), "[[", 1),
    well_position = factor(
      x = well_position,
      levels = unlist(
        lapply(LETTERS[1:16], function(x) paste0(x, 1:24)),
        use.names = TRUE)),
    well_position_in_original_dilution_layout = factor(
      x = cell_type_descriptor,
      levels = unlist(
        lapply(LETTERS[1:16], function(x) paste0(x, 1:24)),
        use.names = TRUE)),
    cell_type_descriptor = NULL,
    sequencing_run = "NN179") %>%
  arrange(plate_number, well_position)

# Construct NN183 sample sheet -------------------------------------------------

file_nn183 <- here(
  "data",
  "sample_sheets",
  "C084_NN183_primer layout_Feb2020.xlsx")

# NOTE: Header row is split across 2 lines, which I combine into 1 before
#       reading in the rest of the spreadsheet.
header_row <- read_excel(
  path = file_nn183,
  sheet = "Sample & Index",
  skip = 2,
  n_max = 1)

header_row <- paste0(colnames(header_row), header_row[1, ])
header_row <- gsub("^\\.\\.\\.[0-9]+", "", header_row)
sample_sheet_nn183 <- read_excel(
  path = file_nn183,
  sheet = "Sample & Index",
  skip = 4,
  col_names = header_row,
  # NOTE: Setting the max guess_max value avoids problems with incorrectly
  #       guessed columns
  #       (https://github.com/tidyverse/readxl/issues/414#issuecomment-352437730)
  guess_max = 1048576)
# Tidy up names
sample_sheet_nn183 <- clean_names(sample_sheet_nn183)
# Remove empty rows/columns.
sample_sheet_nn183 <- remove_empty(sample_sheet_nn183)

# Filter out those empty wells.
sample_sheet_nn183 <- sample_sheet_nn183 %>%
  filter(sample_name != "empty")

# Some final tidying.
sample_sheet_nn183 <- sample_sheet_nn183 %>%
  mutate(
    # NOTE: Some wacky whitespace characters in the plate number
    plate_number = gsub(
      "[\\h\\v]",
      "_",
      trimws(plate_number, whitespace = "[\\h\\v]"),
      perl = TRUE),
    # NOTE: There are some wonky well_positions (e.g., 'I19=A1') that need to
    #       be fixed (these occur because it means well I19 with primer A1,
    #       in SCORE's terminology. I've asked for this to be avoided going
    #       forward.).
    well_position = gsub(" ", "", well_position),
    well_position = sapply(strsplit(well_position, "="), "[[", 1),
    well_position = factor(
      x = well_position,
      levels = unlist(
        lapply(LETTERS[1:16], function(x) paste0(x, 1:24)),
        use.names = TRUE)),
    well_position_in_original_dilution_layout = factor(
      x = cell_type_descriptor_from_original_dilution_layout,
      levels = unlist(
        lapply(LETTERS[1:16], function(x) paste0(x, 1:24)),
        use.names = TRUE)),
    cell_type_descriptor_from_original_dilution_layout = NULL,
    # NOTE: There are some wonky RPIs
    illumina_index_index_number_separate_index_read = trimws(
      illumina_index_index_number_separate_index_read),
    sequencing_run = "NN183") %>%
  arrange(plate_number, well_position)

# Construct sample metadata ----------------------------------------------------

metadata_xls <- here(
  "data", "sample_sheets", "ZIPT SAMPLE LIST WITH PLATE LOCATIONS.xlsx")

baseline_tbl <- read_excel(metadata_xls, range = "B3:K100")
colnames(baseline_tbl) <- paste0(
  colnames(baseline_tbl),
  ifelse(is.na(baseline_tbl[1, ]), "", baseline_tbl[1, ]))
colnames(baseline_tbl) <- gsub("^\\.\\.\\.[0-9]+", "", colnames(baseline_tbl))
baseline_tbl <- baseline_tbl[-1, ]
baseline_tbl <- clean_names(baseline_tbl)
colnames(baseline_tbl) <- sub("cy_tof", "cytof", colnames(baseline_tbl))
baseline_tbl$timepoint <- factor("baseline", levels = c("baseline", "endline"))

endline_tbl <- read_excel(metadata_xls, range = "O3:X104")
colnames(endline_tbl) <- paste0(
  colnames(endline_tbl),
  ifelse(is.na(endline_tbl[1, ]), "", endline_tbl[1, ]))
colnames(endline_tbl) <- gsub("^\\.\\.\\.[0-9]+", "", colnames(endline_tbl))
endline_tbl <- endline_tbl[-1, ]
endline_tbl <- clean_names(endline_tbl)
colnames(endline_tbl) <- sub("cy_tof", "cytof", colnames(endline_tbl))
endline_tbl$timepoint <- factor("endline", levels = c("baseline", "endline"))

metadata_tbl <- rbind(baseline_tbl, endline_tbl) %>%
  mutate(
    collection_date = as.Date(collection_date, tryFormats = c("%d.%m.%Y")),
    sample_serum = ifelse(sample_serum == "√", TRUE, FALSE),
    rna = ifelse(rna == "√", TRUE, FALSE),
    cytof_1 = ifelse(cytof_1 == "√", TRUE, FALSE),
    cytof_2 = ifelse(cytof_2 == "√", TRUE, FALSE))

# Construct final sample sheet -------------------------------------------------

sample_sheet <- rbind(sample_sheet_nn179, sample_sheet_nn183) %>%
  inner_join(
    metadata_tbl,
    by = c("well_position_in_original_dilution_layout" = "plate_location")) %>%
  mutate(
    well_position = factor(
      x = well_position,
      levels = unlist(
        lapply(LETTERS[1:16], function(x) paste0(x, 1:24)),
        use.names = TRUE)),
    well_position_in_original_dilution_layout = factor(
      x = well_position_in_original_dilution_layout,
      levels = unlist(
        lapply(LETTERS[1:16], function(x) paste0(x, 1:24)),
        use.names = TRUE)),
    rowname = paste0(plate_number, "_", well_position)) %>%
  arrange(sequencing_run, plate_number, well_position) %>%
  tibble::column_to_rownames("rowname") %>%
  DataFrame(., check.names = FALSE)

# NOTE: Check that there aren't any malformed well_positions (e.g., 'I19=A1').
stopifnot(!anyNA(sample_sheet$well_position))

# Key variables ----------------------------------------------------------------

plates <- unique(sample_sheet$plate_number)
names(plates) <- plates
sequencing_runs <- tapply(
  sample_sheet$sequencing_run,
  sample_sheet$plate_number,
  unique)
outdir <- here("data", "SCEs")
dir.create(outdir, recursive = TRUE)
extdir <- here("extdata", sequencing_runs, "scPipe", plates)
names(extdir) <- plates
sapply(extdir, dir.create, recursive = TRUE)
# NOTE: Only using first 7 nt of barcode.
read_structure <- get_read_str("CEL-Seq2")
read_structure$bl2 <- 7
# NOTE: Must be an element of biomaRt::listDatasets(), e.g.,
#       biomaRt::listDatasets(biomaRt::useEnsembl("ensembl"))[["dataset"]]
organism <- "hsapiens_gene_ensembl"
# NOTE: Must be an element of biomaRt::listAttributes(), e.g.,
#       biomaRt::listAttributes(biomaRt::useEnsembl("ensembl", organism))[["name"]]
gene_id_type <- "ensembl_gene_id"

# Input files ------------------------------------------------------------------

# FASTQ files
r1_fq <- c(
  grep(
    pattern = "Undetermined",
    x = list.files(
      path = here("extdata", "NN179"),
      full.names = TRUE,
      pattern = glob2rx("*R1.fastq.gz")),
    invert = TRUE,
    value = TRUE),
  grep(
    pattern = "Undetermined",
    x = list.files(
      path = here("extdata", "NN183", "merged"),
      full.names = TRUE,
      pattern = glob2rx("*R1.fastq.gz")),
    invert = TRUE,
    value = TRUE))

r2_fq <- gsub("R1", "R2", r1_fq)
stopifnot(all(file.exists(r2_fq)))

tx_fq <- file.path(extdir, paste0(plates, ".R2.fastq.gz"))
names(tx_fq) <- plates
barcode_fq <- gsub("R2", "R1", tx_fq)

# Concatenate FASTQ files at the plate-level.
mclapply(plates, function(plate) {
  message(plate)
  sequencing_run <- unique(
    sample_sheet[sample_sheet$plate_number == plate, "sequencing_run"])
  if (sequencing_run == "NN179") {
    rpi <- sub(
      " ",
      "-",
      unique(
        sample_sheet[sample_sheet$plate_number == plate,
                     "illumina_index_index_number_separate_index_read"]))
    cmd <- paste0(
      "cat ",
      grep(rpi, r1_fq, value = TRUE),
      " > ",
      barcode_fq[[plate]],
      "\n",
      "cat ",
      grep(rpi, r2_fq, value = TRUE),
      " > ",
      tx_fq[[plate]])
  } else if (sequencing_run == "NN183") {
    rpi <- paste0(
      trimws(
        sub(
          "RPI",
          "",
          unique(
            sample_sheet[sample_sheet$plate_number == plate,
                         "illumina_index_index_number_separate_index_read"]),
          perl = TRUE),
        whitespace = "[\\h\\v]"),
      "_S")
    cmd <- paste0(
      "cat ",
      grep(sequencing_run, grep(rpi, r1_fq, value = TRUE), value = TRUE),
      " > ",
      barcode_fq[[plate]],
      "\n",
      "cat ",
      grep(sequencing_run, grep(rpi, r2_fq, value = TRUE), value = TRUE),
      " > ",
      tx_fq[[plate]])
  }
  system(cmd)
})

# Genome index
genome_index <- here("extdata", "GRCh38.p12", "GRCh38_with_ERCC")

# Genome annotation(s)
annofn <- c(
  here("extdata", "GRCh38.p12", "gencode.v28.primary_assembly.annotation.gff3"),
  system.file("extdata", "ERCC92_anno.gff3", package = "scPipe"))

# Cell barcodes
bc_anno <- file.path(extdir, paste0(plates, ".barcode_annotation.csv"))
names(bc_anno) <- plates

for (plate in plates) {
  message(plate)
  tmp <- sample_sheet[sample_sheet$plate_number == plate, ]
  barcode_df <- data.frame(
    cell_id = row.names(tmp),
    # NOTE: For some reason the primer name and sequence columns have been
    #       reversed in this sample sheet.
    # NOTE: Only using first 7 nt of barcode.
    barcode = strtrim(tmp$c_rt1_primer_name, 7),
    stringsAsFactors = FALSE)
  stopifnot(!anyDuplicated(barcode_df$barcode))
  write.csv(
    x = barcode_df,
    file = bc_anno[[plate]],
    quote = FALSE,
    row.names = FALSE)
}

# Output files -----------------------------------------------------------------

combined_fq <- file.path(extdir, gsub("R[12]", "combined", basename(tx_fq)))
names(combined_fq) <- names(tx_fq)
subread_bam <- gsub("fastq.gz", "subread.bam", combined_fq, fixed = TRUE)
exon_bam <- gsub("subread", "exon", subread_bam)

# FASTQ reformatting -----------------------------------------------------------

filter_settings <- list(rmlow = TRUE, rmN = FALSE, minq = 20, numbq = 2)
# NOTE: Have to loop over files because sc_trim_barcode() is not vectorised.
mclapply(seq_along(tx_fq), function(i) {
  message(combined_fq[i])
  sc_trim_barcode(
    outfq = combined_fq[i],
    r1 = tx_fq[i],
    r2 = barcode_fq[i],
    read_structure = read_structure,
    filter_settings = filter_settings)
})

# Aligning reads to a reference genome -----------------------------------------

Rsubread::align(
  index = genome_index,
  readfile1 = combined_fq,
  output_file = subread_bam,
  nthreads = 12)

# Assigning reads to annotated exons -------------------------------------------

bam_tags <- list(am = "YE", ge = "GE", bc = "BC", mb = "OX")
bc_len <- read_structure$bl1 + read_structure$bl2
barcode_vector <- ""
UMI_len <- read_structure$ul
stnd <- TRUE
fix_chr <- FALSE
mclapply(seq_along(subread_bam), function(i) {
  message(i)
  sc_exon_mapping(
    inbam = subread_bam[i],
    outbam = exon_bam[i],
    annofn = annofn,
    bam_tags = bam_tags,
    bc_len = bc_len,
    barcode_vector = barcode_vector,
    UMI_len = UMI_len,
    stnd = stnd,
    fix_chr = fix_chr)
})

# De-multiplexing data ---------------------------------------------------------

max_mis <- 1
has_UMI <- TRUE
mito <- "chrM"
mclapply(seq_along(exon_bam), function(i) {
  message(i)
  sc_demultiplex(
    inbam = exon_bam[i],
    outdir = extdir[i],
    bc_anno = bc_anno[i],
    max_mis = max_mis,
    bam_tags = bam_tags,
    mito = mito,
    has_UMI = has_UMI)
})

# Gene counting deduped data ---------------------------------------------------

UMI_cor <- 1
gene_fl <- FALSE
mclapply(seq_along(bc_anno), function(i) {
  message(i)
  sc_gene_counting(
    outdir = extdir[i],
    bc_anno = bc_anno[i],
    UMI_cor = UMI_cor,
    gene_fl = gene_fl)
})

# Create and save deduped SingleCellExperiment ---------------------------------

list_of_sce <- lapply(plates, function(plate) {
  create_sce_by_dir(
    datadir = extdir[[plate]],
    organism = organism,
    gene_id_type = gene_id_type,
    pheno_data = sample_sheet[sample_sheet$plate_number == plate, ],
    # NOTE: Create the report separately for more fine-grained control.
    report = FALSE)
})
sce <- Reduce(function(x, y) .combine(x, y, rowData_by = NULL), list_of_sce)
sce$UMI_deduped <- TRUE
assay(sce, withDimnames = FALSE) <- as(
  assay(sce, withDimnames = FALSE),
  "dgCMatrix")
sce <- splitAltExps(sce, ifelse(isSpike(sce), "ERCC", "Endogenous"))
sce <- clearSpikes(sce)

saveRDS(
  sce,
  file.path(outdir, "C084_Hayman_Pasricha.UMI_deduped.scPipe.SCE.rds"),
  compress = "xz")

# Create QC report of deduped data ---------------------------------------------

library(readr)
library(plotly)
library(DT)
library(scater)
library(scran)
library(Rtsne)
# NOTE: Needs a fix for https://github.com/LuyiTian/scPipe/issues/100.
dir.create(here("output", "scPipe"), recursive = TRUE)
# NOTE: Tends to crap itself if using mclapply().
lapply(plates, function(plate) {
  try(create_report(
    sample_name = plate,
    outdir = extdir[[plate]],
    r1 = tx_fq[[plate]],
    r2 = barcode_fq[[plate]],
    outfq = combined_fq[[plate]],
    read_structure = read_structure,
    filter_settings = filter_settings,
    align_bam = subread_bam[[plate]],
    genome_index = genome_index,
    map_bam = exon_bam[[plate]],
    exon_anno = annofn,
    stnd = stnd,
    fix_chr = fix_chr,
    barcode_anno = bc_anno[[plate]],
    max_mis = max_mis,
    UMI_cor = UMI_cor,
    gene_fl = gene_fl,
    organism = organism,
    gene_id_type = gene_id_type))

  # NOTE: Workaround bug in create_report() and stop output after 'Data summary'
  #       section.
  tmp <- readLines(file.path(extdir[[plate]], "report.Rmd"))
  tmp <- c(tmp[1:161], "knitr::knit_exit()", tmp[162:length(tmp)])
  writeLines(tmp, file.path(extdir[[plate]], "report.Rmd"))
  knitr::wrap_rmd(
    file = file.path(extdir[[plate]], "report.Rmd"),
    width = 120,
    backup = NULL)
  rmarkdown::render(
    input = file.path(extdir[[plate]], "report.Rmd"),
    output_file = file.path(extdir[[plate]], "report.html"),
    knit_root_dir = ".")

  # NOTE: Copy the QC report to the repository.
  file.copy(
    from = file.path(extdir[[plate]], "report.nb.html"),
    to = here(
      "output",
      "scPipe",
      paste0(plate, ".scPipe_QC_report.nb.html")),
    overwrite = TRUE)
})

# Gene counting non-deduped data -----------------------------------------------

sapply(paste0(extdir, "_no_dedup"), dir.create)
lapply(names(bc_anno), function(n) {
  message(n)
  geneCountingNoUMIDedup(
    outdir = extdir[[n]],
    bc_anno = bc_anno[[n]])
})

# Create and save non-deduped SingleCellExperiment -----------------------------

list_of_no_dedup_sce <- lapply(plates, function(plate) {
  x <- data.table::fread(
    file.path(paste0(extdir[[plate]], "_no_dedup"), "gene_count.csv"))
  counts <- as.matrix(x[, -1])
  sce <- SingleCellExperiment(
    list(counts = counts),
    colData = sample_sheet[colnames(counts), ])
  rownames(sce) <- x$gene_id
  sce
})
no_dedup_sce <- Reduce(
  function(x, y) .combine(x, y, rowData_by = NULL),
  list_of_no_dedup_sce)
no_dedup_sce$UMI_deduped <- FALSE
colnames(no_dedup_sce) <- paste0(colnames(sce), ".not_UMI_deduped")
assay(no_dedup_sce, withDimnames = FALSE) <- unname(
  as(assay(no_dedup_sce, withDimnames = FALSE), "dgCMatrix"))
no_dedup_sce <- splitAltExps(
  no_dedup_sce,
  ifelse(grepl("^ERCC", rownames(no_dedup_sce)), "ERCC", "Endogenous"))
no_dedup_sce <- clearSpikes(no_dedup_sce)
# NOTE: Order genes and samples as in `sce`.
no_dedup_sce <- no_dedup_sce[rownames(sce),
                             paste0(colnames(sce), ".not_UMI_deduped")]
saveRDS(
  no_dedup_sce,
  file.path(outdir, "C084_Hayman_Pasricha.not_UMI_deduped.scPipe.SCE.rds"),
  compress = "xz")
