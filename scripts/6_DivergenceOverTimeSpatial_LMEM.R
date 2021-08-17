### OC SK 2021
# This script runs the linear mixture model for divergence from ancestral strain
# be that tp1, consensus M, or consensus C
# per 1000 bp's
# input fasta per patient with tp's and aligned ancestor.

### ref
# https://ourcodingclub.github.io/tutorials/mixed-models/


### options
analysis_outdir = "analysis/6_DivergenceOverTimeLMEM/"
indir = "data/tree/"
patients = c("15664","16207","22763","22828","26892","28545","29447","47939")
vl_file = "data/patient_vl.csv" # we need thsis to match tp's to dates
bps = 1000 # how wide to look?
###

#----- functions
library(lme4)
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


# ----- generate data
starts = c(1,1000,2000,3000,4000,5000,6000,7000,8000)
for(start in starts){ # spatial
  end = start + bps
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
    
    
    # decide which sample to calculate diverence from
    # 0 sequence is an infereed ancestor
    # 1 sequence is the first tp
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
    # remove initial tp
    df = df[df$months != 0,]
    
    df2 = rbind(df2,df)
    
    df2 = df2[df2$tp1_div < 0.4,] # remove spurious
    

    
  }

}



# clean 

# ----- Linear Mixed Effects Modelling

# clean up the explanatory variable so centred and scaled
df2$norm_months= scale(df2$months, center = TRUE, scale = TRUE)


mixed.lmer <- lmer(tp1_div ~ norm_months + (1|patient), data = df2)
summary(mixed.lmer)

plot(mixed.lmer)  # looks alright, no patterns evident


qqnorm(resid(mixed.lmer))
qqline(resid(mixed.lmer))  # points fall nicely onto the line - good!







# ----- plot the LMEM per patient
mixed.ranslope <- lmer(tp1_div ~ norm_months + (1 + norm_months|patient), data = df2)

### plot
ggplot(df2, aes(x = norm_months, y = tp1_div, colour = as.factor(patient))) +
  geom_point(alpha = 0.5) +
  theme_classic() +
  geom_line(data = cbind(df2, pred = predict(mixed.ranslope)), aes(y = pred), size = 1)   # adding predicted line from mixed model 
  #theme(legend.position = "none",
  #      panel.spacing = unit(2, "lines"))  # adding space between panels




# ----- plot the LMEM overall

library(ggeffects)

# Extract the prediction data frame
pred.mm <- ggpredict(mixed.lmer, terms = c("norm_months"))  # this gives overall predictions for the model

# Plot the predictions 

(ggplot(pred.mm) + 
    geom_line(aes(x = x, y = predicted)) +          # slope
    geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                fill = "lightgrey", alpha = 0.5) +  # error band
    geom_point(data = df2,                      # adding the raw data (scaled values)
               aes(x = norm_months, y = tp1_div, colour = as.factor(patient))) + 
    labs(x = "norm_months", y = "divergence from tp1 sample", 
         title = "Virus may be diverging over time") + 
    theme_minimal()
)



# extract stats
library(stargazer)
stargazer(mixed.lmer, type = "text",
          digits = 3,
          star.cutoffs = c(0.05, 0.01, 0.001),
          digit.separator = "")
