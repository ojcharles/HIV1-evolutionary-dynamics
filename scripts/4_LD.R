library(ggrepel)
library(wesanderson)
library(ggplot2)


fasta2ld = function(alignment, outname){
  plink = "plink1.9"
  #msa -> vcf
  command = paste0("wsl snp-sites ", alignment, " -v -o ", alignment, ".vcf")
  system(command)

  #vcf -> plink
  command = paste0("wsl ", plink, " --vcf ",alignment,".vcf --double-id --out data/ld/", outname)
  system(command)

  # calculate LD
  command = paste0("wsl ", plink, " --bfile data/ld/",outname,"  --r2 dprime --inter-chr --ld-window-r2 0.0 --out data/ld/", outname)
  system(command)


  # C:\Oscar\OneDrive\UCL\21-2_hiv_tasp\data\LD>python WeightedLD.py --input 16207.fasta --min-acgt 0.99 --min-variability 0.1 --r2-threshold 0.0 >
  
}

fasta2ld(alignment, "one")



# this script generates the LD figures for the paper

# these were altered a lot

# the actual LD processing from fasta to Ld matrix was done with weighted LD to allow control of acgt fraction, i ran these weighted
# we can ignore that in the paper as it has very little effect here.



ld_all = data.frame()


# -----------------------plot all R2 values
patients = c("15664", "16207", "22763")
for(patient in patients){
  dat.ld = read.table(paste0("data/LD/",patient,".wld"), header = T, sep = "")
  dat.ld$range = dat.ld$site_b  - dat.ld$site_a
  
  png(paste0("analysis/4_ld/",patient, "_ld_hist.png"))
  hist(dat.ld$r2,breaks = 100,main = paste(patient, "- LD distribution, minor allele freq 0.1, min acgt 0.99"))
  dev.off()
  pal <- wes_palette("Zissou1", 100, type = "continuous")
  
  
  g = ggplot(data = dat.ld[dat.ld$r2 >= 0 ,]) +
    geom_point(aes(x=site_a, y =site_b, colour = r2),size = 0.1) +
    scale_colour_gradientn(colours = pal) +
    theme(
      # get rid of panel grids
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      legend.direction = "horizontal",
      legend.position='bottom',
    ) + theme_classic()
  
  
  ggsave(plot = g,filename = paste0("analysis/4_ld/",patient, "_ld_spatial.png"), device = "png")
  
  
  
  ########## LD DECAY
  
  
  # mean values
  bases = 250
  th = data.frame()
  for(i in seq(1,max(dat.ld$site_a) - (bases + 100), bases)){ # ignore end bases
    r2_vals = dat.ld[dat.ld$range > i & dat.ld$range < (i + bases) & dat.ld$site_b < 7000 ,]$r2 # what are the ld values between range x and y
    #r2_vals = dat.ld[dat.ld$range > i & dat.ld$range < (i + bases) ,]$r2
    th = rbind(th,data.frame(dist = i, r2_mean = mean(r2_vals)))
  }
  th$patient = patient
  ld_all = rbind(ld_all, th)

}



g = ggplot(ld_all,aes(x = dist, y = r2_mean, colour = patient)) +
  geom_point(alpha = 0.2) +
  geom_smooth(span = 0.5) +
  theme_classic()

ggsave("analysis/4_ld/3_pats_ld_decay.png", g,"png")




