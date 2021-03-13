# plot all patients
# all variants in grey
# drm variants in colour with legend
# no facet



#-------- load all functions into environment
funfiles = list.files("R",full.names = T)
sapply(funfiles, function(x) source(x))




#-------- vcf -> vcf.csv
files_vcf = list.files("1", pattern = "*.vcf$",recursive = T,full.names = T)
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




#========= plots
dirs = list.dirs("1")[-1]
patients = str_split(dirs, "/", simplify = T)[,2]

patients_which = c("15664", "16207", "22763")
patients_which = c(3,4,8)

for(i in patients_which){
  i = 3
  print(paste(i, "------", patients[i]))
  
  files = list.files(dirs[i], pattern = "*.vcf.csv", full.names = T)
  file = files[1]
  
  dat = read.csv(file,fileEncoding="UTF-8-BOM", stringsAsFactors = F)[,-1]
  dat$freq = readr::parse_number(dat$freq)
  
  dat.all = data.frame(timepoint = 1, mutation = dat$change, freq = dat$freq, consequence = dat$CONSEQUENCE)
  
  for(k in 2:length(files)){
    
    file = files[k]
    dat = read.csv(file,fileEncoding="UTF-8-BOM", stringsAsFactors = F)[,-1]
    dat$freq = readr::parse_number(dat$freq)
    
    t.dat = data.frame(timepoint = k, mutation = dat$change, freq = dat$freq, consequence = dat$CONSEQUENCE)
    
    
    dat.all = rbind(dat.all, t.dat)
    
  }
  
  dat.all$mutation = as.character(dat.all$mutation)
  dat.all$consequence = as.character(dat.all$consequence)
  
  ## clean
  dat.plot = dat.all[dat.all$consequence == "synonymous",]
  for(mut in unique(dat.plot$mutation)){
    if(min(dat.plot[dat.plot$mutation ==mut,]$freq) >= 95.0){ # remove if all freq above 95%
      dat.plot = dat.plot[dat.plot$mutation != mut,]
    }
    else if(max(dat.plot[dat.plot$mutation ==mut,]$freq) <= 5.0){ # remove if all freq below 5%
      dat.plot = dat.plot[dat.plot$mutation != mut,]
    }
    else if(length(dat.plot[dat.plot$mutation ==mut,]$freq) == 1 )# remove if only 1 freq
      dat.plot = dat.plot[dat.plot$mutation != mut,]
  }
  
  dat.plot$gene = str_split(dat.plot$mutation, "_",simplify = T)[,1]
  dat.plot$mut = str_split(dat.plot$mutation, "_",simplify = T)[,2]
  dat.plot = dat.plot[dat.plot$gene %in% c("rt", "integrase", "protease"),]
  
 
  
  
  ## plots
  
  # ggplotly()
  
  
  ### update - grey out non res plots
  dat.drm = read_csv("1/major_drms.csv")
  dat.drm$change = paste(dat.drm$gene, dat.drm$mutation, sep = "_")
  dat.drm.plot = dat.all[dat.all$consequence == "nonsynonymous",]
  dat.drm.plot$gene = str_split(dat.drm.plot$mutation, "_",simplify = T)[,1]
  dat.drm.plot$mut = str_split(dat.drm.plot$mutation, "_",simplify = T)[,2]
  dat.drm.plot = dat.drm.plot[dat.drm.plot$mutation %in% dat.drm$change,]
  
  ## - fill in empty synonymous
  for(tp in 1:max(dat.plot$timepoint)){
    for(mut in dat.plot$mutation)
      if(nrow(dat.plot[dat.plot$mutation == mut & dat.plot$timepoint == tp,]) == 0)
        dat.plot = rbind(dat.plot, data.frame(timepoint = tp, mutation = mut,  freq = 0,   consequence = "",
                                                      gene = str_split(mut, "_",simplify = T)[,1], mut = str_split(mut, "_",simplify = T)[,2]))
  }
  # remove dupliacte synonymous
  for(mut in dat.plot$mutation){
    if(nrow(dat.plot[dat.plot$mutation == mut,]) > max(dat.plot$timepoint)){ # if mutations has any dupliacates
      dat.plot = dat.plot[dat.plot$mutation != mut,]
    }
  }

  
  
  
  
  ## - fill in empty drm
  if(nrow(dat.drm.plot) > 0){
    for(tp in 1:max(dat.drm.plot$timepoint)){
      for(mut in dat.drm.plot$mutation)
        if(nrow(dat.drm.plot[dat.drm.plot$mutation == mut & dat.drm.plot$timepoint == tp,]) != 1)
          dat.drm.plot = rbind(dat.drm.plot, data.frame(timepoint = tp, mutation = mut,  freq = 0,   consequence = "",
                                         gene = str_split(mut, "_",simplify = T)[,1], mut = str_split(mut, "_",simplify = T)[,2]))
    }
  }
  
  g = ggplot(data = dat.plot, aes(x = timepoint, y = freq)) +
    geom_line(aes(group = mutation), alpha = 0.1) + # syn
    geom_line(data = dat.drm.plot, aes( colour = mutation), alpha = 0.8, size = 1) + # drm
    geom_point(data = dat.drm.plot, aes( colour = mutation), alpha = 0.8, size = 2) +
    #facet_grid( ~ gene) +
    theme_classic() +
    labs(title = str_split(basename(file), "_", simplify = T)[1])
    
  plot(g)
  
  outfile = paste("1/plots/", str_split(basename(file), "_", simplify = T)[1], ".png", sep = "" )
  ggsave(outfile, device = "png",width = 15, height = 10)
  
  
}# dirs

