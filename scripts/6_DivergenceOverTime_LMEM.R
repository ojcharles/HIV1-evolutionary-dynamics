### OC SK 2021
# This script runs the linear mixture model for divergence from ancestral strain
# be that tp1, consensus M, or consensus C
# over the whole genome
# input fasta per patient with tp's and aligned ancestor.

### ref
# https://ourcodingclub.github.io/tutorials/mixed-models/


### options
analysis_outdir = "analysis/6_DivergenceOverTime/"
indir = "data/tree/"
patients = c("15664","16207","22763","22828","26892","28545","29447","47939")
vl_file = "data/patient_vl.csv" # we need thsis to match tp's to dates
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
library(ggeffects)
library(stargazer)
library(lmerTest)
elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
} # takes timepoints, matches to dates in the VL data. then normalises such that we get months from first datapoint


# ----- generate data
ancestries = c("aligned_baseline_" , "aligned_c_", "aligned_m_")
for(ancestry in ancestries){
  
  #starts = c(1,1000,2000,3000,4000,5000,6000,7000,8000)
  starts = 1
  for(start in starts){
    df2 = data.frame() # rowbind into
    for(i in 1:length(patients)){
      patient = patients[i]
      alignment = paste0(indir, ancestry, patient, ".fasta")
      msa = read.dna(alignment,format = "fasta",as.matrix = T) #read msa
      end = length(msa[1,])
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
      
      
      # sequence of index 1 is the inferred ancestor, so we just want the diversity per timepoint from this
      col = grep("1|timepoint",colnames(df)) # get col indexes of interest
      df = df[,col] # subset cols
      colnames(df) = c("diversity", "timepoint")
      
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
      
      df2 = df2[df2$diversity < 0.4,] # remove spurious
      
      
    } # patient
  } # start
  
  
  
  
  
  # ----- per ancestry run this
  # ----- Linear Mixed Effects Modelling
  # is there an association between diversity and time, when we control for the variation in each patient?
  
  # clean up the explanatory variable so centred and scaled
  df2$norm_months= scale(df2$months, center = TRUE, scale = TRUE)
  
  mixed.lmer <- lme4::lmer(diversity ~ months + (1|patient), data = df2)
  summary(mixed.lmer)
  
  plot(mixed.lmer)  # looks alright, no patterns evident
  
  qqnorm(resid(mixed.lmer))
  qqline(resid(mixed.lmer))  # points fall nicely onto the line - good!
  
  
  
  
  # ----- plot the LMEM per patient
  mixed.ranslope <- lme4::lmer(diversity ~ months + (1 + norm_months|patient), data = df2)
  
  ### plot
  g = ggplot(df2, aes(x = months, y = diversity, colour = as.factor(patient))) +
    geom_point(alpha = 0.5) +
    theme_classic() +
    geom_line(data = cbind(df2, pred = predict(mixed.ranslope)), aes(y = pred), size = 1)   # adding predicted line from mixed model 
  
  outfile = paste0(analysis_outdir, "WG_",ancestry, "random-slope-random-intercept-model.png")
  ggsave(outfile,g,"png")
  
  
  
  # ----- plot the LMEM overall
  
  # Extract the prediction data frame
  pred.mm <- ggpredict(mixed.lmer, terms = c("months"))  # this gives overall predictions for the model
  
  # Plot the predictions 
  g = (ggplot(pred.mm) + 
      geom_line(aes(x = x, y = predicted)) +          # slope
      geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                  fill = "lightgrey", alpha = 0.5) +  # error band
      geom_point(data = df2,                      # adding the raw data (scaled values)
                 aes(x = months, y = diversity, colour = as.factor(patient))) + 
      labs(x = "norm_months", y = "divergence from tp1 sample", 
           title = "Virus may be diverging over time") + 
      theme_minimal()
  )
  outfile = paste0(analysis_outdir, "WG_",ancestry, "overall-model-prediction.png")
  ggsave(outfile,g,"png")
  
  
  
  
  
  # ----- extract stats
  # table
  stargazer(mixed.lmer, type = "text",
            digits = 3,
            star.cutoffs = c(0.05, 0.01, 0.001),
            digit.separator = "")
  
  # confidence intervals on the effect of months
  t = confint(mixed.lmer)
  print(t) # if the months 2.5 and 97.5 are both + or both -  then that's fun.
  
  
  
  
} # ancestry

file.remove("del.fasta")

















