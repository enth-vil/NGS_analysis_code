ref <- read.csv("~/Desktop/test_pupil/PANCAN_PDAC_100plex_ref.csv")
library(tidyverse)
bed <- ref %>%
  select(Chr, Absolute.Start, Absolute.End, Gene, strand) 
library("GenomicRanges")
library(rtracklayer) 
bedgr <- GRanges(bed[,1], IRanges(bed[,2], bed[,3]), bed$strand)
export(bedgr, "filter_pos.bed", "bed")








bed <- ref %>%
  separate(ID, into = c("name_chrom", "chromStart_End"), sep = ":", convert = TRUE) %>%
  separate(chromStart_End , into = c("chromStart" , "chromEnd"), sep = "-", convert = TRUE) %>%
  separate(name_chrom , into = c("name" , "chrom"), sep = "_", convert = TRUE) %>%
  select(chrom, chromStart, chromEnd, name, strand)
library("GenomicRanges")
library(rtracklayer)
bedgr <- GRanges(bed[,1], IRanges(bed[,2], bed[,3]), bed$strand)
export(bedgr, "filter_pos.bed", "bed")

