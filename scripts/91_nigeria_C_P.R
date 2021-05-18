# this is a script to cluster each gene fasta file & calculate KNN accuracy as a score of how well the gene separates the earupean & african sequences

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

# genome scan
# how many clusters along genome
files = list.files(path = "data/nigeria_C_P/",pattern = "*.fasta$",full.names = T)

for(i in 1:length(files)){

  
  
  # read and churn data
  
  file = files[i]
  seq <- read.dna(file = file, format = "fasta",as.matrix = T)
  labels = stringr::str_split(labels(seq),"_",simplify = T)
  t.seq = seq[2:nrow(seq),]
  
  # filtering
  {
  #emove = c("12796", "13003", "17339", "19937", "22097", "23271", "29219", "31254", "34253")
  #remove = "28841"
  #remove = which(grepl(remove, labels[,1]))
  #keep = setdiff(1:length(labels[,1]), remove)
  #t.seq = t.seq[keep,]
  }
  
  t.dist <- dist.dna(t.seq, model = "TN93", pairwise.deletion = T, as.matrix = T)
  # get contribtion of axis http://r-sig-ecology.471788.n2.nabble.com/variance-explained-for-PCoA-axes-td6239213.html
  t.mds <- cmdscale(t.dist, k = 2, eig = T)
  t.axis_contr = t.mds$GOF
  t.mds = t.mds$points
  t.mds = data.frame(sample = rownames(t.mds), MDS_D1 = t.mds[,1], MDS_D2 = t.mds[,2])
  
  
  
  # optimal clusters - using kmeans & a consensus method
  ######################## NBClust provides 30 indices for determining the optimal number of clusters, a but like randomforest
  # so rather than BIC or AIC or AIC3, this is a process that sums ovar all methods.
  nb <- NbClust(t.mds[,2:3], diss=NULL, distance = "euclidean",
                min.nc=2, max.nc=9, method = "kmeans",
                index = "all", alphaBeale = 0.1)
  
  cluster = as.data.frame(nb$Best.partition)
  cluster$sample = rownames(cluster)
  colnames(cluster) = c("cluster", "sample")
  cluster$cluster = as.factor(cluster$cluster)
  df = merge(x = t.mds,y = cluster,by = "sample")
  
  # what to colour by
  df$group = str_split(df$sample, "_", simplify = T)[,2]
  df$group[df$group == "2C"] <- "C"
  df$group[df$group == "C"] <- "CSF"
  df$group[df$group == "P"] <- "Plasma"
  
  
  write.csv(df, paste0(file, ".csv"))
  
  
  
  base = 6 # set the height of your figure (and font)
  expand = 2 # font size increase (arbitrarily set at 2 for the moment)
  
  g = ggplot(df,aes(x = MDS_D1, y = MDS_D2, colour = group, label = sample),alpha = 0.5) +
    geom_point(size = 4, alpha = 0.8) +
    theme_pubr() +
    labs(x = "", y = "",title = paste0("Patient: 00", i)) +
    theme(legend.position="right") +
    guides(colour=guide_legend(title="Sample type")) +
    theme(axis.text.x = element_text(angle = 90, hjust = 0, size = base * expand ),
          axis.text.y = element_text(size = base * expand )) +
    scale_color_manual(values = c("Plasma" = "#F8766D",
                                  "CSF"="#00BFC4"))
  
  plot(g)
  ggsave(plot = g, paste0(file, ".png"), device = "png",dpi = 1000)
  p = plotly::ggplotly(g)
  htmlwidgets::saveWidget(p,file.path(normalizePath(dirname(paste0(file, ".html"))),basename(paste0(file, ".html"))))
}
