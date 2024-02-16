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

# ### process
readr::read_tsv(input_file) %>%
  dplyr::select(lane, sample_name, barcode) %>%
  dplyr::arrange(sample_name) %>%
  tidyr::nest(-lane) %>%
  purrr::walk2(.x = .$data, .y = .$lane, .f = ~ readr::write_tsv(.x, paste0(.y, ".txt"), col_names = FALSE))

lanes <- readr::read_tsv(input_file) %>%
  .$lane %>%
  unique

for(i in 1:length(lanes)){
  file_buff <- readr::read_tsv(paste0(lanes[i], ".txt"), col_names = F)
  file.create(paste0(lanes[i], ".fasta"))
  for(x in 1:nrow(file_buff)){
    id<-paste0(">", file_buff$X1[x])
    barcode<-paste0("^", file_buff$X2[x])
    readr::write_lines(paste0(lanes[i], ".fasta"), x=id, append = T)
    readr::write_lines(paste0(lanes[i], ".fasta"), x=barcode, append = T)
  }
}
