# script to read vcf files, and return the csv files.
# this is to check we converted them right

library(hivdrg)

#-------- vcf -> vcf.csv
files_vcf = list.files("data/vcf-2/", pattern = "*.vcf$",recursive = T,full.names = T)
# only intereted in what happens to dodgy files

file = "data/vcf-2/15664/15664_2.vcf"

for(file in files_vcf){
  outfile = paste(file, ".csv", sep = "")
  if(file.exists(outfile) == FALSE){
    # handle data
    dat1 = read_input(file, global)
    ### ANNOTATE VARIANTS
    dat2 <- annotate_variants(f.dat = dat1, global)
    #edit
    #---
    write.csv(dat2, outfile)
  }
}


# 184.8