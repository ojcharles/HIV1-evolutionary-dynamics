# this script takes the patient separated fast files, as used in the phylogeny and generates the pairwise distances from tp1 to all other tps.
# Then plots this divergence per month with smoothing



### options
analysis_outdir = "analysis/6_DivergenceOverTime/"
indir = "data/tree/"
patients = c("15664", "16207","22763")
vl_file = "data/patient_vl.csv"
###

#----- functions
library(ggplot2)
library(stringr)
library(ggnewscale)
library(RColorBrewer)
library(gridExtra)
library(dplyr)
library(ape)
elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
} # takes timepoints, matches to dates in the VL data. then normalises such that we get months from first datapoint



df2 = data.frame() # rowbind into
for(i in 1:length(patients)){
  patient = patients[i]
  
  msa = read.FASTA(paste0(indir, patient, ".fasta")) # read alignment
  dist = dist.dna(msa,model = "TN93", pairwise.deletion = TRUE, as.matrix = T) # generate pairwise distances
  df = data.frame(dist)
  # clean labels and just leave tp
  tp = str_split(rownames(df),"_",simplify = T)[,2]
  rownames(df) = tp
  tp = str_split(colnames(df),"_",simplify = T)[,2]
  colnames(df) = tp
  df$timepoint = tp
  
  # filter as we just want distance from tp1
  col = grep("1|timepoint",colnames(df))
  df = df[,col]
  colnames(df) = c("tp1_div", "timepoint")
  
  # we now have a vector we are sure represents all sequences divergence from tp1
  

#---------------------------- read in tp - months data
{
  vl_dat = read.csv("data/patient_vl.csv")
  vl_dat = vl_dat[ vl_dat$patient == patients[i] , ]
  ###### careful the line bwelow is a large loop
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
}
  
  
# timepoints -> months
df = merge(df,vl_dat, by = "timepoint")

df2 = rbind(df2,df)
}



#------------- simple R default smoothing

g = ggplot(df2,aes(x = months, y = tp1_div, colour = as.factor(patient))) +
  geom_point() +
  geom_smooth(se = F) +
  theme_classic() +
  labs(subtitle = "linear regression of TN93 pairwise distannces of all timepoints against the patient tp0 sequence as proxy for founder sequence") +
  guides(colour=guide_legend(title="Patient"))
#g
ggsave(filename = paste0(analysis_outdir, "divergence_smooth.png"), device = "png",plot = g)

ggsave(filename = paste0(analysis_outdir, "divergence_smooth.tiff"), device = "tiff",plot = g)





#------------- linear model

g = ggplot(df2,aes(x = months, y = tp1_div, colour = as.factor(patient))) +
  geom_point() +
  geom_smooth(method = "lm",se = F) +
  theme_classic() +
  labs(subtitle = "linear regression of TN93 pairwise distannces of all timepoints against the patient tp0 sequence as proxy for founder sequence") +
  guides(colour=guide_legend(title="Patient"))
#g
ggsave(filename = paste0(analysis_outdir, "divergence_regression.png"), device = "png",plot = g)

ggsave(filename = paste0(analysis_outdir, "divergence_regression.tiff"), device = "tiff",plot = g)

