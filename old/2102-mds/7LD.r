# Linkage Disequilibium

# using plink 1.9
#https://www.cog-genomics.org/plink/1.9/
# ld stuff https://www.cog-genomics.org/plink/1.9/ld
plink = "/home/oscar/tools/plink"
  
# # plink 2 https://www.cog-genomics.org/plink/2.0/
# plink = "/home/oscar/tools/plink2"


# ----- generate vcf from alignment & generate homologous plot
command = paste("snp-sites", alignment, "-v")
system(command)

# move to folder
command = paste("mv *.vcf 7-mutinf/")
system(command)



# ----- vcf to LD

command = paste0(plink, " --vcf 7-mutinf/",basename(alignment),".vcf --double-id --out 7-mutinf/plink") 
system(command)


command = paste(plink, "--bfile 7-mutinf/plink --r2 --inter-chr --mac 5 --out 7-mutinf/plink-ld") # minor allelle of 5 freq
system(command)


# waheyy it works!
# plot
dat.ld = read.table("7-mutinf/plink-ld.ld", header = T, sep = "")

g = ggplot(data = dat.ld) +
  geom_tile(aes(x=BP_A, y =BP_B, colour = R2)) +
  scale_color_gradient2(low = "white", high = "red") + 
  theme(
    # get rid of panel grids
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change plot and panel background
    plot.background=element_rect(fill = "gray"),
    panel.background = element_rect(fill = 'black'),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = "black", color = NA),
    legend.key = element_rect(color = "gray", fill = "black"),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white")
  )
ggsave(g,filename = "7-mutinf/plink-ld.png", device = "png",width = 20, height = 20)


# ----- vcf to LD  biallellic only

command = paste0(plink, " --vcf 7-mutinf/",basename(alignment),".vcf --double-id --biallelic-only --out 7-mutinf/plinkba") # re-running with biallellic only
system(command)


command = paste(plink, "--bfile 7-mutinf/plinkba --r2 --inter-chr --mac 5 --out 7-mutinf/plinkba-ld")
system(command)


# waheyy it works!
# plot
dat.ld = read.table("7-mutinf/plinkba-ld.ld", header = T, sep = "")

g = ggplot(data = dat.ld) +
  geom_tile(aes(x=BP_A, y =BP_B, colour = R2)) +
  scale_color_gradient2(low = "white", high = "red") + 
  theme(
    # get rid of panel grids
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    # Change plot and panel background
    plot.background=element_rect(fill = "gray"),
    panel.background = element_rect(fill = 'black'),
    legend.direction = "horizontal",
    legend.background = element_rect(fill = "black", color = NA),
    legend.key = element_rect(color = "gray", fill = "black"),
    legend.title = element_text(color = "white"),
    legend.text = element_text(color = "white")
  )
ggsave(g,filename = "7-mutinf/plinkba-ld.png", device = "png",width = 20, height = 20)

