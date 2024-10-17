#!/usr/bin/env Rscript

################################################################################
# combine processed counts by summing counts by sample name across lanes
# part of CRISPR / shRNA screen pre-processing pipeline
# 
# Jesse J. Lipp
# Institute of Molecular Pathology (IMP), Vienna, Austria
# 2017/09/20
################################################################################

### command line parameters
args         <- commandArgs(trailingOnly = TRUE)
library_file <- args[1]
count_files  <- args[2:length(args)]

### functions
`%>%` <- dplyr::`%>%`

read_featurecounts <- function(path) {
  readr::read_tsv(path, skip = 1) %>%
    dplyr::select(-Chr, -Start, -End, -Strand, -Length) %>%
    dplyr::rename(id = Geneid)
}

### combine counts
library <- readr::read_tsv(library_file) %>%
  dplyr::select(id, group)

names(count_files) <- stringr::str_replace(basename(count_files), ".txt", "")
pattern <- paste(c(paste0(names(count_files), "_"), ".sam"), collapse = "|")

lapply(count_files, read_featurecounts) %>%
  purrr::map(tidyr::gather, sample_name, count, -id) %>%
  dplyr::bind_rows() %>%
  dplyr::mutate(sample_name = stringr::str_replace_all(sample_name, pattern, "")) %>%
  #remove stagger length information from sample name
  dplyr::mutate(sample_name = stringr::str_split(sample_name, "#") %>% lapply("[[", 2) %>% unlist %>%
                  stringr::str_split("STAGGERLENGTH") %>% lapply("[[", 1) %>% unlist) %>%
  dplyr::group_by(id, sample_name) %>%
  dplyr::summarize(count = sum(count)) %>%
  dplyr::ungroup() %>%
  tidyr::spread(sample_name, count) %>%
  dplyr::inner_join(library, by = "id") %>%
  dplyr::select(id, group, dplyr::everything()) %>%
  readr::format_tsv() %>%
  cat
