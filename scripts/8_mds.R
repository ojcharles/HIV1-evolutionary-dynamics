#Haplotype plot

library(ape)              # distance analysis & sikmple DNA object
library(Biostrings)       # Provides DNAString, DNAStringSet, etc
library(GenomicRanges)    # Provides GRanges, etc
library(GenomicFeatures)  # provides Txdb
library(seqinr)           # why this?
library(ggplot2)
library(stringr)
library(NbClust)
library(class)
library(ggrepel)
library(mclust)
library(ggpubr)
library(reshape2)
library(vegan)
set.seed(123)


# ----- inputs
file = "data/all_patients_haplos.fasta"
file = "data/LD/n=8consensus_msa.fasta"
mds_method = "non-metric" # metric or non-metric



# read in sequence data to distance matrix
t.seq <- read.dna(file = file, format = "fasta",as.matrix = T)
t.dist <- dist.dna(t.seq, model = "TN93", pairwise.deletion = T, as.matrix = T)
labels = stringr::str_split(labels(seq),"_",simplify = T)

# MDS
if(mds_method == "metric"){
  # normal metric like pca
  t.mds <- cmdscale(t.dist, k = 2, eig = T)
  t.axis_contr = t.mds$GOF
  t.mds = t.mds$points
  t.mds = data.frame(sample = rownames(t.mds), MDS_D1 = t.mds[,1], MDS_D2 = t.mds[,2])
}else{
  # non metric
  t.mds = metaMDS(t.dist, k=2, try = 1000, trymax = 5000, maxit = 5000, plot = T)
  t.mds = t.mds$points
  t.mds = data.frame(sample = rownames(t.mds), MDS_D1 = t.mds[,1], MDS_D2 = t.mds[,2])
  stressplot(t.mds) # if the points are fairly around the red then we are happy
}



### this code assigns each sequence to a cluster based on Kmeans
# optimal clusters - using kmeans & a consensus method
nb <- NbClust(t.mds[,2:3], diss=NULL, distance = "euclidean",
              min.nc=2, max.nc=9, method = "kmeans",
              index = "all", alphaBeale = 0.1)
cluster = as.data.frame(nb$Best.partition)
cluster$sample = rownames(cluster)
colnames(cluster) = c("cluster", "sample")
cluster$cluster = as.factor(cluster$cluster)
df = merge(x = t.mds,y = cluster,by = "sample")
df$patient = str_split(df$sample, "_", simplify = T)[,1]


# title = "figure n shows the multi dimensional scaling of every haplotypes whole genome pairwise distances using a TN93 substitution model
# into the two dimensions of maximum variance. An ensemble model using 30 indices has been used to determine the number of clusters and their assignment.
# ref: https://www.jstatsoft.org/article/view/v061i06 "

write.csv(df, paste0("analysis/8_mds/",basename(file), ".clusters.csv"))

man_col_pal = c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7",
                "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD",
                "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D",
                "#8A7C64", "#599861")

base = 6 # set the height of your figure (and font)
expand = 2 # font size increase (arbitrarily set at 2 for the moment)

g = ggplot(df,aes(x = MDS_D1, y = MDS_D2, colour = patient, label = sample),alpha = 0.5) +
  geom_point(size = 4, alpha = 0.8) +
  theme_pubr() +
  labs(x = "", y = "", title = paste("consensus for all patients -", mds_method, " MDS method" )) +
  theme(legend.position="right") +
  guides(colour=guide_legend(title="Patient")) +
  scale_colour_manual(values = man_col_pal) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0, size = base * expand ),
        axis.text.y = element_text(size = base * expand ))


ggsave(plot = g, paste0("analysis/8_mds/",basename(file), ".png"), device = "png",dpi = 1000)
ggsave(plot = g, paste0("analysis/8_mds/",basename(file), ".tiff"), device = "tiff",dpi = 300,width = 8, height = 8)
# p = plotly::ggplotly(g)
# htmlwidgets::saveWidget(p,file.path(normalizePath(dirname(paste0(file, ".html"))),basename(paste0(file, ".html"))))


