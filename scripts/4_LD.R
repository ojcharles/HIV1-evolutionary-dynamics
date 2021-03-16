setwd("/home/oscar/steve-hiv")


library(ggrepel)

plink = "/home/oscar/tools/plink"
alignment = "LD/consensus_filt_LD.fasta"

# ----- generate vcf from alignment & generate homologous plot
command = paste("snp-sites", alignment, "-v")
system(command)

# move to folder
command = paste("mv *.vcf LD/")
system(command)


# ----- vcf to LD

command = paste0(plink, " --vcf LD/",basename(alignment),".vcf --double-id --out LD/plink") 
system(command)


command = paste(plink, "--bfile LD/plink --r2 --inter-chr --mac 5 --out LD/plink-ld") # minor allelle of 5 freq
system(command)

# -----------------------plot all R2 values
dat.ld = read.table("LD/plink-ld.ld", header = T, sep = "")

g = ggplot(data = dat.ld) +
  geom_tile(aes(x=BP_A, y =BP_B, colour = R2)) +
  scale_color_gradient2(low = "white", high = "red") + 
  theme(
    # get rid of panel grids
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.direction = "horizontal",
    legend.position='bottom',
  ) + theme_classic()
ggsave(g,filename = "LD/ld.png", device = "png",width = 20, height = 20)





# -----------------------plot R2 above 0.5
dat.ld = read.table("LD/plink-ld.ld", header = T, sep = "")
dat.ld = dat.ld[dat.ld$R2 >= 0.5,]

g = ggplot(data = dat.ld) +
  geom_tile(aes(x=BP_A, y =BP_B, colour = R2)) +
  scale_color_gradient2(low = "white", high = "red") + 
  theme(
    # get rid of panel grids
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.direction = "horizontal",
    legend.position='bottom',
  ) + theme_classic()

ggsave(g,filename = "LD/ld2.png", device = "png",width = 20, height = 20)
