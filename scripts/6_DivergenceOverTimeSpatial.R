# this script takes the patient separated fast files, as used in the phylogeny and generates the pairwise distances from tp1 to all other tps.
# Then plots this divergence per month with smoothing
# now does this per 500 base pairs to look at spatial behaviour, does this do fun stuff at the highly LD sites?



### options
analysis_outdir = "analysis/6_DivergenceOverTime/"
indir = "data/tree/"
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
library(ape)
library(ggpubr)
elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
} # takes timepoints, matches to dates in the VL data. then normalises such that we get months from first datapoint


#starts = 1000
starts = c(1,1000,2000,3000,4000,5000,6000,7000,8000)
for(start in starts){ # spatial
  end = start + 1000
  df2 = data.frame() # rowbind into
  for(i in 1:length(patients)){
    patient = patients[i]
    msa = read.dna(paste0(indir, patient, ".fasta"),format = "fasta",as.matrix = T) #read msa
    msa = msa[,start:end]
    write.FASTA(msa, "del.fasta")
    msa = read.FASTA("del.fasta") # read submsa
    dist = dist.dna(msa,model = "TN93", pairwise.deletion = TRUE, as.matrix = T) # generate pairwise distances
    df = data.frame(dist)
    # clean labels and just leave tp
    tp = str_split(rownames(df),"_",simplify = T)[,2]
    rownames(df) = tp
    tp = str_split(colnames(df),"_",simplify = T)[,2]
    colnames(df) = tp
    df$timepoint = tp
    
    # filter as we just want distance from ancestral reconstructed subtype C (same in all patients)
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
  
  df2 = df2[df2$tp1_div <= 0.4,] # remove outliers
  #------------- per patient plot of data + linear model
  
  my.formula <- y ~ x
  
  g = ggplot(df2,aes(x = months, y = tp1_div, colour = as.factor(patient))) +
    geom_point() +
    geom_smooth(method = "lm", se = F, formula = my.formula) +
    theme_classic() + theme(text = element_text(size=26)) +
    scale_x_continuous(breaks = c(0,5,10,17,21,27,31,33)) + 
    xlab("Months since switch to 2nd-Line Regimen") +
    ylab("Diversity") +
    labs(subtitle = paste("all patient", start, ":", end)) + # remove?
    guides(colour=guide_legend(title="Patient"))+ # remove?
    stat_cor() 
  
  g
  ggsave(filename = paste0(analysis_outdir, "divergence_regression",start,".png"), device = "png",plot = g)
  
  #ggsave(filename = paste0(analysis_outdir, "divergence_regression",patient,".tiff"), device = "tiff",plot = g)
}


#------------- simple R default smoothing

#g = ggplot(df2,aes(x = months, y = tp1_div, colour = as.factor(patient))) +
# geom_point() +
#geom_smooth(se = F) +
#theme_classic() +
#labs(subtitle = "linear regression of TN93 pairwise distannces of all timepoints against the patient tp0 sequence as proxy for founder sequence") +
#guides(colour=guide_legend(title="Patient"))
#g
#ggsave(filename = paste0(analysis_outdir, "divergence_smooth.png"), device = "png",plot = g)

#ggsave(filename = paste0(analysis_outdir, "divergence_smooth.tiff"), device = "tiff",plot = g)


#------------- linear model

g = ggplot(df2,aes(x = months, y = tp1_div, colour = as.factor(patient))) +
  geom_point() +
  geom_smooth(method = "lm", se = T,) +
  theme_classic(text = element_text(size=26)) +
  scale_x_continuous(breaks = c(0,5,10,17,21,27,31,33)) + 
  xlab("Months since switch to 2nd-Line Regimen") +
  ylab("Diversity") +
  guides(colour=guide_legend(title="Patient"))

#g
ggsave(filename = paste0(analysis_outdir, "divergence_regression.png"), device = "png",plot = g)

ggsave(filename = paste0(analysis_outdir, "divergence_regression.tiff"), device = "tiff",plot = g)

