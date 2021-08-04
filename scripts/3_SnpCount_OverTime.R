# March 2021 - data for limitations and such
# per patient - per timepoint generates the number of Synonymous and nonynonymous mutations
# uses filters as in script 2 for variant dynamics

### options
analysis_outdir = "analysis/3_SnpCount_OverTime/"
indir = "data/vcf-mar21-trust"
patients = c("15664","16207","22763","22828","26892","28545","29447","47939")
vl_file = "data/patient_vl.csv"
###



#----- functions
library(ggplot2)
library(stringr)
library(ggnewscale)
library(RColorBrewer)
library(gridExtra)
library(dplyr)

elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
} # takes timepoints, matches to dates in the VL data. then normalises such that we get months from first datapoint

get_filtered_data_per_consequence <- function(dat.all, consequence, vl_file){
  dat.plot = dat.all[dat.all$consequence == consequence,c(1:4, 7,6)]
  for(mut in unique(dat.plot$mutation)){
    if(length(unique(dat.plot[dat.plot$mutation ==mut,]$loc)) > 1) { # remove if two separate loci -> AA change are detected take the one with highest max value
      t = dat.plot[dat.plot$mutation ==mut,]
      t.keep = t$loc[which.max(t$freq)] # keep
      #print(t)
      t.throw = unique(t$loc)[unique(t$loc)!= t.keep]
      dat.plot = dat.plot[dat.plot$loc != t.throw,]
    } 
    
    if(min(dat.plot[dat.plot$mutation ==mut,]$freq) >= 95.0){ # remove if all freq above 95%  
      dat.plot = dat.plot[dat.plot$mutation != mut,]
    }
    else if(max(dat.plot[dat.plot$mutation ==mut,]$freq) <= 5.0){ # remove if all freq below 5%
      dat.plot = dat.plot[dat.plot$mutation != mut,]
    }
    else if(length(dat.plot[dat.plot$mutation ==mut,]$freq) == 1 ){# remove if only 1 freq
      dat.plot = dat.plot[dat.plot$mutation != mut,]
    }
  }
  
  dat.plot$gene = str_split(dat.plot$mutation, "_",simplify = T)[,1]
  dat.plot$mut = str_split(dat.plot$mutation, "_",simplify = T)[,2]
  dat.plot = dat.plot[dat.plot$gene %in% c("rt", "integrase", "protease", "gag", "env"),]
  
  
  
  #-------------------------------------------------------- patient specific filters, ie remove any syn for patient a at tp 2 between 1k and 2k as no coverage
  if(patients[i] == 31254){
    dat.plot = dat.plot[dat.plot$timepoint != 4,]
  }
  if(patients[i] == 28545){
    dat.plot = dat.plot[dat.plot$timepoint != 2,]
  }
  # updates march 2021
  if(patients[i] == 15664){
    # remove, all ENV synonymous changes at positions 5849 - 9758 at timepoint 4 (month 21)
    dat.plot = dat.plot[!(dat.plot$timepoint == 4 & dat.plot$loc %in% 5849:9758),]
  }
  if(patients[i] == 22763){
    # remove, all ENV synonymous changes at positions 5849 - 9758 at timepoint 4 (month 21)
    dat.plot = dat.plot[!(dat.plot$timepoint == 2 & dat.plot$loc %in% 789:1907 ),]
  }
  
  
  
  
  
  #---------------------------- now add column that links timepoint to date
  vl_dat = read.csv("data/patient_vl.csv")
  vl_dat = vl_dat[ vl_dat$patient == patients[i] , ]
  ###### careful the line bwelow is a large loop
  if(nrow(vl_dat) == 0){ # if there is no vl data
    outfile = paste("1/plots1.1/", str_split(basename(file), "_", simplify = T)[1], ".png", sep = "" )
    ggsave(plot = g, filename = outfile, device = "png",width = 15, height = 10)
    next
  }
  
  
  vl_dat$timepoint = 1:nrow(vl_dat)
  vl_dat$sample_date = as.Date(vl_dat$sample_date, "%d/%m/%Y")
  #for paper
  for(j in 1:nrow(vl_dat)){
    if(j > 1){
      vl_dat$months[j] = elapsed_months(vl_dat$sample_date[j], vl_dat$sample_date[1])
    }else{vl_dat$months[j] = 0}
  }
  # forced bodge
  if(patient == "16207"){
    vl_dat[vl_dat$months == 4,]$months <- 3
  }
  dat.plot = merge(dat.plot, vl_dat, by = "timepoint")
  
  
  return(dat.plot)
} # this mimics exactly the filtering processes in script 1

collapse_per_tp <- function(df){
  out = df %>% group_by(months) %>% tally
  return(out)
} # collapses table of varfreqs per months -> count of freqs per month. 




df.s = data.frame() ; df.ns = data.frame() # setting up tables for all pat plot later 
iter = 1:length(patients)
for(i in iter){
  patient = patients[i]
  print(paste(i, "------", patients[i]))
  
  # read in vcf.csv files, keep only those related to patient
  files = list.files(indir, pattern = "*.vcf.csv", full.names = T)
  files = files[grep(patient, files)]
  file = files[1]
  
  # load first file
  dat = read.csv(file,fileEncoding="UTF-8-BOM", stringsAsFactors = F)[,-1]
  dat$freq = readr::parse_number(dat$freq)
  dat.all = data.frame(timepoint = 1, mutation = dat$change, freq = dat$freq, consequence = dat$CONSEQUENCE, dna_ref = dat$REFCODON, dna_var = dat$VARCODON, loc = dat$start)
  
  # append all timepoints varaints files
  for(k in 2:length(files)){
    file = files[k]
    dat = read.csv(file,fileEncoding="UTF-8-BOM", stringsAsFactors = F)[,-1]
    dat$freq = readr::parse_number(dat$freq)
    t.dat = data.frame(timepoint = k, mutation = dat$change, freq = dat$freq, consequence = dat$CONSEQUENCE, dna_ref = dat$REFCODON, dna_var = dat$VARCODON, loc = dat$start)
    dat.all = rbind(dat.all, t.dat)
  }
  
  # clean    
  dat.all$mutation = as.character(dat.all$mutation)
  dat.all$consequence = as.character(dat.all$consequence)
  dat.all$gene = str_split(dat.all$mutation, "_",simplify = T)[,1]
  dat.all$mut = str_split(dat.all$mutation, "_",simplify = T)[,2]
  dat.all = dat.all[dat.all$gene %in% c("rt", "integrase", "protease","gag", "env"),]

  cons = "synonymous"
  d = get_filtered_data_per_consequence(dat.all, cons, vl_file)
  d = collapse_per_tp(d)
  write.csv(d,paste0(analysis_outdir,patient, cons, ".csv"))
  df.s = rbind(df.s,data.frame(patient = patient, d))
  
  cons = "nonsynonymous"
  d = get_filtered_data_per_consequence(dat.all, cons, vl_file)
  d = collapse_per_tp(d)
  write.csv(d,paste0(analysis_outdir,patient, cons, ".csv"))
  df.ns = rbind(df.ns, data.frame(patient = patient, d))
}


for(patient in unique(df.s$patient)){
  tp0 = max(df.s[df.s$patient == patient & df.s$months == 0,]$n)
  df.s[df.s$patient == patient,3] = df.s[df.s$patient == patient,3] - tp0
}




#---------------- plot all together
ggplot(df.ns[!df.s$patient == 22763,],aes(x = months, y = n)) +
  geom_point(aes( colour = patient)) +
  geom_smooth(method = "lm") +
  ggpubr::stat_cor()






