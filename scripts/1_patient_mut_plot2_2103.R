# March 2021 - plots for three aptients Steve wants for the preprint
# variant dynamics
# resistance mutations
# timepoint for VL

### options
analysis_outdir = "analysis/2103-patientvcf2data/"
indir = "data/patientvcf/"
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

# colour palette
mycolours <- c(
  "dodgerblue2", "#E31A1C", # red
  "green4",
  "#6A3D9A", # purple
  "#FF7F00", # orange
  "black",
  "maroon", "orchid1", "blue1",
  "darkturquoise", "green1", "#ff00ff", "#6495ed",
  "darkorange4"
)
mycolours = RColorBrewer::brewer.pal(12,"Paired") # dark2



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
  
  
  
  #-------------------------------------------------------- Variant plot
  
  #palette = "Manual"
  g = ggplot(data = dat.plot, aes(x = months, y = freq)) +
    #geom_line(aes(group = mutation, colour = gene), alpha = 0.2, size = 1) + # syn
    geom_line(aes(group = mutation, colour = `Genome Position`), alpha = 0.2) + # syn
    #scale_colour_gradientn(colours = terrain.colors(10)) +
    #scale_color_gradient(low = "#0080ff", high = "black") +
    scale_color_gradient2(midpoint = 5000, low = "red", mid = "green",
                          high = "blue", space = "Lab" ) +
    #scale_color_gradientn(colours = rainbow(3)) +
    #scale_color_gradient() +
    new_scale_colour() + # required as now two colour axis
    geom_line(data = dat.drm.plot, aes( colour = mutation), alpha = 0.8, size = 1) + # drm
    geom_point(data = dat.drm.plot, aes( colour = mutation), alpha = 0.8, size = 2, position=position_jitter(w=0, h=1)) +
    #facet_grid( ~ gene) +
    theme_classic() +
    # labs(subtitle = paste(str_split(basename(file), "_", simplify = T)[1], "jittered", palette)) +
    scale_y_continuous(breaks = c(0,20,40,60,80,100)) +
    scale_x_continuous(labels = (unique(vl_dat$months)),breaks=c(unique(vl_dat$months))) +
    #scale_x_continuous(breaks = (unique(vl_dat$months)), labels=c(0,4,5,9,12,17)) +
    #scale_x_continuous(breaks = (unique(vl_dat$months)), labels=c(0,4,18,19,25,28)) +
    labs(x = "Time since baseline (Months)", y = "Mutationa Frequency (%)") +
    #scale_color_brewer(palette=palette) 
    # manual colour scale
    #guides(colour=guide_legend(title="Major DRMs"))
    scale_color_manual(values=sample(mycolours), name = "Major DRMs")
  
  # headers - automatic from steves notes
  dat.header = data.frame(patient = c("22828", "22763", "16207","15664","29219","31254"), header = c("Patient 22828, Female, 28yrs",
                                                                                                     "Patient 22763, Female, 42yrs",
                                                                                                     "Patient 16207, Female, 35yrs",
                                                                                                     "Patient 15664, Male, 51yrs",
                                                                                                     "Patient 29219, something something",
                                                                                                     "Patient 31254, something something"))
  if(patients[i] %in% dat.header$patient){
    g = g + labs(subtitle = dat.header$header[which(dat.header$patient == patients[i])])}
  
  plot(g)
  
  
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
  
  g2
  
  
  g_all = grid.arrange(g2,g, ncol = 1, heights = c(1, 4))
  
  
  
  outfile = paste0(analysis_outdir, patient, ".png" )
  ggsave(plot = g_all, filename = outfile, device = "png",width = 15, height = 10)
  outfile = paste0(analysis_outdir, patient, ".tiff" )
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






