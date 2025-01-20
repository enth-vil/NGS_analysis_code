ref <- read.csv("~/Desktop/test_pupil/PANCAN_PDAC_100plex_ref.csv")
library(tidyverse)
bed <- ref %>%
  select(Chr, Absolute.Start, Absolute.End, Gene, strand) 
bed$Chr <-  gsub("chr","",as.character(bed$Chr))
library("GenomicRanges")
library(rtracklayer) 
bedgr <- GRanges(bed[,1], IRanges(bed[,2], bed[,3]), bed$strand)
export(bedgr, "filter_pos.bed", "bed")
