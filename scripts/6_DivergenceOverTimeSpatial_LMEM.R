### OC SK 2021
# This script runs the linear mixture model for divergence from ancestral strain
# be that tp1, consensus M, or consensus C
# for 1000 bp chunks
# input fasta per patient with tp's and aligned ancestor.

### ref
# https://ourcodingclub.github.io/tutorials/mixed-models/
# https://cran.r-project.org/web/packages/lme4/vignettes/lmer.pdf


### options
analysis_outdir = "analysis/6_DivergenceOverTime/"
indir = "data/tree/"
patients = c("15664","16207","22763","22828","26892","28545","29447","47939")
vl_file = "data/patient_vl.csv" # we need thsis to match tp's to dates
bps = 1000 # how wide to look?
CILevel # confidence Interval with FDR. so we want 95%, but doing 9 tests. so a bonferroni level would be
###

#-------------------- functions
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
library(Mediana) # adjust CI for FDR
elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
} # takes timepoints, matches to dates in the VL data. then normalises such that we get months from first datapoint




# -------------------- per founder, per genomic region, calculate: lmem & p value for effect of months
ancestries = c("aligned_baseline_" , "aligned_c_", "aligned_m_")
for(ancestry in ancestries){
  
  # -----  per genomic section
  starts = c(1,1000,2000,3000,4000,5000,6000,7000,8000)
  #starts = 1
  dfp = data.frame() # p values for effect of months per chunk
  for(start in starts){
    # ----- generate data
    df2 = data.frame() # x and y values
    for(i in 1:length(patients)){
      patient = patients[i]
      alignment = paste0(indir, ancestry, patient, ".fasta")
      msa = read.dna(alignment,format = "fasta",as.matrix = T) #read msa
      end = start + bps
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
      if(length(grep("0", colnames(df))) > 0){ # if there is a tp0 patient
        keep_tp1 = T
        col = grep("0|timepoint",colnames(df)) # get col indexes of interest
      }else{
        keep_tp1 = F
        col = grep("1|timepoint",colnames(df)) # get col indexes of interest
      }
      col = grep("1|timepoint",colnames(df)) # get col indexes of interest
      df = df[,col] # subset cols
      colnames(df) = c("diversity", "timepoint")
      
      # we now have a vector we are sure represents all sequences divergence from tp1
      
      
      # ----- read in tp - months data
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
      if(keep_tp1 == F){df = df[df$months != 0,]}
      
      df2 = rbind(df2,df)
      
      df2 = df2[df2$diversity < 0.4,] # remove spurious
      
    } # patient
    
    # remove NA, these occure where its all just insertions
    df2 = df2[complete.cases(df2),]
  
    
    
    
    
    # ---------- Stats with Linear Mixed Effects Modelling
    # we are testing for is the slope (effect) of months on diversity significantly different from 0?
    # + or -
    
    # clean up the explanatory variable so centred and scaled
    # we actually ran this with normalised and erbatim months and the same thing
    df2$norm_months= scale(df2$months, center = TRUE, scale = TRUE)
    
    # ----- generate the model with an overall trend 
    # lmerTest is a pull of lmer4 with p value stats (compared to 0) using a Satterthwaite t-test
    mixed.lmer <- lmerTest::lmer(diversity ~ months + (1|patient), data = df2)
    summary(mixed.lmer)
    # dev
    #plot(mixed.lmer)  # looks alright, no patterns evident
    #qqnorm(resid(mixed.lmer))
    #qqline(resid(mixed.lmer))  # points fall nicely onto the line - good!
    
    
    
    
    
    # ----- plot the LMEM per patient
    # give it a random slope per patient
    mixed.ranslope <- lmerTest::lmer(diversity ~ months + (1 + months|patient), data = df2)
    
    ### plot
    g1 = ggplot(df2, aes(x = months, y = diversity, colour = as.factor(patient))) +
      geom_point(alpha = 0.5) +
      theme_classic() +
      geom_line(data = cbind(df2, pred = predict(mixed.ranslope)), aes(y = pred), size = 1)   # adding predicted line from mixed model 
    #geom_line(data = cbind(df2, pred = predict(mixed.lmer)), aes(y = pred), size = 1)   # yep so the above is a random slope
    # I beleive the mixed.randslope model after running ggplotly
    outfile1 = paste0(analysis_outdir, start,"_",ancestry, "random-slope-random-intercept-model.png")
    
    
    
    
    # ----- plot the LMEM overall
    
    # Extract the prediction data frame
    pred.mm <- ggpredict(mixed.lmer, terms = c("months"))  # this gives overall predictions for the model
    
    # Plot the predictions 
    g2 = (ggplot(pred.mm) + 
            geom_line(aes(x = x, y = predicted)) +          # slope
            geom_ribbon(aes(x = x, ymin = predicted - std.error, ymax = predicted + std.error), 
                        fill = "lightgrey", alpha = 0.5) +  # error band
            geom_point(data = df2,                       # adding the raw data (scaled values)
                       aes(x = months, y = diversity, colour = as.factor(patient))) + 
            labs(x = "months", y = "Divergence from 1st sample") + 
            theme_minimal(base_size = 14))
    
    outfile2 = paste0(analysis_outdir, start, "_",ancestry, "overall-model-prediction.png")
    
    
    # plot them all
    ggsave(outfile1,g1,"png")
    ggsave(outfile2,g2,"png")
    
    
    
    
    
    
    # ---------- KEY STATISTICS
    
    # what is the p that the effect of months is non-0?
    p = anova(mixed.lmer)$`Pr(>F)`
    

    # now record the p per chunk
    dfp = rbind(dfp, data.frame(start = start, 
                                p = p))
    
    # remove the model and outputs to ensure they are recalculated
    rm("mixed.lmer", "p")
  } # start
  # we now have for the founder, per genome chunk the p value
  
  # adjust p values
  dfp$p
  dfp$p_bh = p.adjust(dfp$p,method = "BH")
  write.csv(dfp, paste0(analysis_outdir, ancestry, "_Pvals_BH_adjust.csv"))
  
  
  
  
  
  
  
  # --- now remove any plots that have insignificant p
  # i could have stored them in a list but eh this works
  
  for(start in starts){
    is_signif = dfp[dfp$start == start,3] < 0.05
    if(is_signif == FALSE){
      # remove the plots
      file.remove(paste0(analysis_outdir, start,"_",ancestry, "random-slope-random-intercept-model.png"))
      file.remove(paste0(analysis_outdir, start, "_",ancestry, "overall-model-prediction.png"))
    }
  }
      


  
} # ancestry

file.remove("del.fasta")










# ---------- code graveyard
# # ----- extract stats
# # table
# p = stargazer(mixed.lmer, type = "text",
#               digits = 3,
#               star.cutoffs = c(0.05, 0.05 / 9),
#               digit.separator = "")
# p.fdr = grepl("\\*\\*", p[7]) # if signif grabs it
# 
# # grab the significance
# # confidence intervals on the effect of months
# out2 = confint(mixed.lmer)

# now loop through the significant


