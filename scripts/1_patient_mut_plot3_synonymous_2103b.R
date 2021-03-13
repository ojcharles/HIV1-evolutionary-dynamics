# March 2021 - plots for three aptients Steve wants for the preprint
# variant dynamics
# resistance mutations
# timepoint for VL

# V3 plot nonsynonymous mutations in teh background rather than synonymous
# This version plots synonymous changes in the background
### options
analysis_outdir = "analysis/2103-patientvcf2plot/"
indir = "data/vcf-mar21-trust"
patients = c("15664", "16207","22763")

###



#----- functions
library(ggplot2)
library(stringr)
library(ggnewscale)
library(RColorBrewer)
library(gridExtra)
elapsed_months <- function(end_date, start_date) {
  ed <- as.POSIXlt(end_date)
  sd <- as.POSIXlt(start_date)
  12 * (ed$year - sd$year) + (ed$mon - sd$mon)
}

#colour palette
#mycolours = RColorBrewer::brewer.pal(7,"Dark2") # dark2
mycolours = c(
  #"#2f4f4f", #darkslategray
  "#008000", #green
  "#C1292E", #red # or set to just "red
  "deeppink1", #safety yellow
  "#0000ff", #blue
  "#ff8c00", #darkorange
  "#6A3D9A", #purple
  "#020100") #black




# # a manual discrete df for lines in the final plots
# mycolours_forced = data.frame(mut = c("rt_K103N", "rt_M184V", "rt_P225H", "rt_K238T", "rt_V75I","integrase_G163R", 
#                                       "rt_K103N",                         "protease_l23I", "integrase_E157Q", "protease_M46I",
#                                       "rt_K103N", "rt_M184V", "rt_P225H", "rt_D67N", "rt_K219Q", "rt_K65R", "rt_K70R"),
#                               cols = c("#008000", "#C1292E","#0000ff", "deeppink1", "#ff8c00", "#6A3D9A",
#                                        "#008000", "deeppink1", "#ff8c00", "#6A3D9A",
#                                        "#008000", "#C1292E","#0000ff", "deeppink1", "#ff8c00", "#6A3D9A", "#020100"))
# mycolours_forced = unique(mycolours_forced)
cols = c("rt_K103N" = "#008000",
         "rt_M184V" = "#C1292E",
         "rt_P225H" = "#0000ff",
         "rt_K238T" = "deeppink1",
         "rt_V75I"  = "#ff8c00",
         "integrase_G163R" = "#6A3D9A",
         "protease_L23I" = "deeppink1",
         "integrase_E157Q" = "#ff8c00",
         "protease_M46I" = "#6A3D9A",
         "rt_D67N" = "deeppink1",
         "rt_K219Q" = "#ff8c00",
         "rt_K65R" = "#6A3D9A",
         "rt_K70R" = "#020100")




#========= plots
# dat.duplicate = data.frame(patient = "1", timpoint = 1, mutation = "", freq = 0.0, consensus = "", dna_ref = "", dna_var = "", gene = "", mut = "",stringsAsFactors = F)[-1,]
# dat.drm.all = data.frame(patient = "1", mutation = "")[-1,]
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
  
  
  
  
  #-------------------------------------------------------- prepare synonymous mutations for ploting
  # remove duplicate synonymous - as multiple mutations can result in an AA change. and we only care about AA change. also only in res genes
  dat.plot = dat.all[dat.all$consequence == "synonymous",c(1:4, 7,6)]
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
    # else if(length(dat.plot[dat.plot$timepoint == 1 & dat.plot$mutation == mut,3]) ==1){ # remove mutations if they are more than 50% at tp1
    #   if(dat.plot[dat.plot$timepoint == 1 & dat.plot$mutation == mut,3] > 50){ 
    #     dat.plot = dat.plot[dat.plot$mutation != mut,]
    #   }
    # }
    # ## mar21 synonymous - if a synonymous meets criteria then keep else throw
    # t9 = dat.plot[dat.plot$mutation == mut & dat.plot$consequence == "synonymous",]
    # remove = max(t9$freq) - min(t9$freq)       # remove a mutations if the difference of max and min freq is less than 340%, ie. boring
    # if(remove <40){dat.plot = dat.plot[dat.plot$mutation != mut,]}
  
  }

  # dat.plot is the synonymous dataset for plotting
  dat.plot$gene = str_split(dat.plot$mutation, "_",simplify = T)[,1]
  dat.plot$mut = str_split(dat.plot$mutation, "_",simplify = T)[,2]
  dat.plot = dat.plot[dat.plot$gene %in% c("rt", "integrase", "protease", "gag", "env"),]
  
  #print(length(unique(dat.plot$mutation)))
  
  #-------------------------------------------------------- prepare Drug Resistant Mutations dfor plotting
  # read DRM table
  dat.drm = read.table("data/resmuts2.tsv",sep = "\t",header = T)
  dat.drm$change = paste(dat.drm$gene, dat.drm$mutation, sep = "_")
  dat.drm.plot = dat.all[dat.all$consequence == "nonsynonymous",1:4]
  dat.drm.plot$gene = str_split(dat.drm.plot$mutation, "_",simplify = T)[,1]
  dat.drm.plot$mut = str_split(dat.drm.plot$mutation, "_",simplify = T)[,2]
  dat.drm.plot = dat.drm.plot[dat.drm.plot$mutation %in% dat.drm$change,] # if any non-synnymous changes are in the table, only keep those else remove
  # remove any mutations that have a max % of less than 2%
  for(mut in unique(dat.drm.plot$mutation)){
    t = dat.drm.plot[dat.drm.plot$mutation == mut,]
    if(max(t$freq) <2){
      dat.drm.plot = dat.drm.plot[!dat.drm.plot$mutation == mut,]
    }
  }
  
  ## - fill in empty synonymous - currently may only appear at 1 tp. need to show as 0% if missing value
  for(tp in 1:max(dat.plot$timepoint)){
    for(mut in dat.plot$mutation)
      if(nrow(dat.plot[dat.plot$mutation == mut & dat.plot$timepoint == tp,]) == 0)
        dat.plot = rbind(dat.plot, data.frame(timepoint = tp, mutation = mut,  freq = 0,   consequence = "", loc = max(dat.plot[dat.plot$mutation == mut,5]),
                                              dna_var = dat.plot[dat.plot$mutation == mut,6][1],gene = str_split(mut, "_",simplify = T)[,1], mut = str_split(mut, "_",simplify = T)[,2]
        )) # fudge but works as all the same value
  }
  
  
  
  ## If there are no DRM for a patient then clean up the data so we show no mutations - else itll fall over
  if(nrow(dat.drm.plot) > 0){
    for(tp in 1:max(dat.drm.plot$timepoint)){
      for(mut in dat.drm.plot$mutation)
        if(nrow(dat.drm.plot[dat.drm.plot$mutation == mut & dat.drm.plot$timepoint == tp,]) != 1)
          dat.drm.plot = rbind(dat.drm.plot, data.frame(timepoint = tp, mutation = mut,  freq = 0,   consequence = "",
                                                        gene = str_split(mut, "_",simplify = T)[,1], mut = str_split(mut, "_",simplify = T)[,2]))
    }
  }
  

  
  
  
  #-------------------------------------------------------- patient specific filters, ie remove any syn for patient a at tp 2 between 1k and 2k as no coverage
  if(patients[i] == 31254){
    dat.plot = dat.plot[dat.plot$timepoint != 4,]
    dat.drm.plot = dat.drm.plot[dat.drm.plot$timepoint != 4,]
  }
  if(patients[i] == 28545){
    dat.plot = dat.plot[dat.plot$timepoint != 2,]
    dat.drm.plot = dat.drm.plot[dat.drm.plot$timepoint != 2,]
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
  
  
  
  
  
  
  mid = 5000
  dat.plot$`Genome Position` = dat.plot$loc

  #-------------------------------------------------------- prepare viral load data for plotting
  
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
  
  
  ## edit
  #if(i == 4){vl_dat$months[2] = 3}
  #vl_dat$months[3] = 16
  ## end
  dat.plot = merge(dat.plot, vl_dat, by = "timepoint")
  dat.drm.plot = merge(dat.drm.plot, vl_dat, by = "timepoint")
  
  # march 2021 remove 5 mutations 
  if(patient == "22763"){
    remove = c("integrase_T97A", "protease_M46I", "protease_M46L", "rt_A98G", "rt_E138G")
    dat.drm.plot = dat.drm.plot[!dat.drm.plot$mutation %in% remove, ] # keep
  }
  if(patient == "16207"){
    remove = c("protease_M46L")
    dat.drm.plot = dat.drm.plot[!dat.drm.plot$mutation %in% remove, ] # keep
  }
  
  
  
  #-------------------------------------------------------- Variant plot
  # t = dat.drm.plot
  # for(r in 1:nrow(t)){
  #   t$col[r] = mycolours_forced[mycolours_forced$mut == t$mutation[r],2]
  #   t$line[r] = mycolours_forced[mycolours_forced$mut == t$mutation[r],3]
  # }
  # dat.drm.plot = t
  
  #palette = "Manual"
  g = ggplot(data = dat.plot, aes(x = months, y = freq)) +
    geom_line(aes(group = mutation, colour = `Genome Position`), alpha = 0.2) + # syn
    scale_color_gradient2(midpoint = 5000, low = "red", mid = "green",
                          high = "blue", space = "Lab" ) +
    theme_classic() +
    scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
    scale_x_continuous(labels = (unique(vl_dat$months)),breaks=c(unique(vl_dat$months))) +
    labs(x = "Time since baseline (Months)", y = "Mutationa Frequency (%)") 
  
  # headers - automatic from steves notes
  dat.header = data.frame(patient = c("22828", "22763", "16207","15664","29219","31254"), header = c("Patient 22828, Female, 28yrs",
                                                                                                     "Patient 22763, Female, 42yrs",
                                                                                                     "Patient 16207, Female, 35yrs",
                                                                                                     "Patient 15664, Male, 51yrs",
                                                                                                     "Patient 29219, something something",
                                                                                                     "Patient 31254, something something"))
  if(patients[i] %in% dat.header$patient){
    g = g + labs(subtitle = dat.header$header[which(dat.header$patient == patients[i])])}
  
  #plot(g)
  
  
  if(patients[i] == 31254){
    vl_dat = vl_dat[vl_dat$timepoint != 4,]
  }
  if(patients[i] == 28545){
    vl_dat = vl_dat[vl_dat$timepoint != 2,]
  }
  
  
  
  
  
  #-------------------------------------------------------- Viral Load Plot
  
  g2 = ggplot(data = vl_dat, aes(x = months, y = log10(vl))) +
    geom_line(linetype = "dashed") +
    geom_point() +
    theme_classic() +
    theme(plot.margin = margin(0,4.0,0,0.4, "cm"),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.x = element_blank())+
    labs(y = "Viral load (log copies/ml)") +
    scale_y_continuous(breaks = seq(0,max(log10(vl_dat$vl)) , by = 0.5)) +
    scale_x_continuous(labels = (unique(vl_dat$months)),breaks=c(unique(vl_dat$months))) +
    #scale_x_continuous(breaks = (unique(vl_dat$months)), labels=c(0,4,5,9,12,17)) +
    #scale_x_continuous(breaks = (unique(vl_dat$months)), labels=c(0,4,18,19,25,28)) +
    expand_limits(y = 0)+
    geom_vline(aes(xintercept=as.numeric(unique(vl_dat$months))))
  
  #g2
  
  
  g_all = grid.arrange(g2,g, ncol = 1, heights = c(1, 4))
  
  
  
  outfile = paste0(analysis_outdir, patient, "synonymous.png" )
  ggsave(plot = g_all, filename = outfile, device = "png",width = 15, height = 10)
  outfile = paste0(analysis_outdir, patient, "synonymous.tiff" )
  ggsave(plot = g_all, filename = outfile, device = "tiff",width = 15, height = 10)
  
  #
  #dat.drm.all = rbind(dat.drm.all, data.frame(patient = patients[i], mutation = unique(dat.drm.plot[,2])))

}# dirs

#write.csv(dat.duplicate, "1/dat_duplicate.csv")
#write.csv(dat.drm.all, "1/dat_drm_all.csv")







# ## drm plot all patients
# library(dplyr)
# dat.drm.all.plot = dat.drm.all %>% 
#   count(mutation) %>% 
#   mutate(perc = n / 18)
# dat.drm.all.plot$gene = factor(str_split(dat.drm.all.plot$mutation, "_", simplify = T)[,1])
# 
# ggplot(dat.drm.all.plot, aes(x = mutation, y = perc, fill = gene)) +
#   geom_bar(stat = "identity") +
#   theme(axis.text.x = element_text(angle = 90)) +
#   facet_wrap(~gene)






