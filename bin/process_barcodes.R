#!/usr/bin/env Rscript

################################################################################
# process sample annotation to obtain input demultiplexing
# part of CRISPR / shRNA screen pre-processing pipeline
# 
# Jesse J. Lipp
# Institute of Molecular Pathology (IMP), Vienna, Austria
# 2017/09/20
################################################################################

# barcode file should be tab-separated text file with three columns:
# 1) lane
# 2) sample_name
# 3) barcode

### command line parameters
args       <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]

### functions
`%>%` <- dplyr::`%>%`

### process
readr::read_tsv(input_file) %>%
  dplyr::select(lane, sample_name, barcode) %>%
  dplyr::arrange(sample_name) %>%
  tidyr::nest(-lane) %>%
  purrr::walk2(.x = .$data, .y = .$lane, .f = ~ readr::write_tsv(.x, paste0(.y, ".txt"), col_names = FALSE))
