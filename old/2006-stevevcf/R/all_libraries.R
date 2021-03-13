# load all function
library(ggplot2)
library(dplyr)
library(stringr)
library(IRanges)
library(Biostrings)
library(GenomicRanges) 
library(GenomicFeatures) # if save txdb as sql may not need to load this
library(VariantAnnotation)
library(reshape2)
library(DT)
library(plotly)
library(readr)


#detach(package:Biostrings)

#package variables
global = list()
global$res_table = "db/judb5.csv"
#create unique session folder
global$date <- format(Sys.time(), "%Y-%m-%d")
global$dir = paste("")

# Steve
global$genome = genome="K03455.1"
global$path_gff3_file="K03455.1.gff3"
global$path_fasta_file="K03455.1.fasta"