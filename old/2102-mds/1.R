# this script runs the whole project

library(hmmcluster)
# devtools::install_github("ojcharles/RCircos") # RCircos in cran falls over when you have overlapping gene regions as blocks
library(RCircos)
library(ggplot2)
library(ape)
library(stringr)
library(ggpubr)
library(aricode)

alignment = "allhaplos.fasta"
msa = ape::read.dna(alignment, format = "fasta",as.matrix = T,as.character = T)
msabin = ape::read.dna(alignment, format = "fasta")

#----------------------------- get regions

# remove the initial repeat region from analysis - 6836 and before should be removed


hmmcluster::get_regions_parallel(alignment, ref.file = "AF411967.fasta", run_width = 30,run_cores = 12,ref.pattern = "FULL",run_model = "BIC")

regions = read.csv("out/3-clean/gt-regions.csv")
regions$width = regions$end - regions$start

# check if any region overlap
for(i in 2:max(regions$region)){
  if(regions$start[i] <= regions$end[i-1]){
    print(i)
  }
}
# add gene information
annotate_gene = function(dat, ignore_genes = c("RL8A", "UL89")){
  #get gene data
  gff3 = system.file("ref/NC_006273.2.gff3", package = "hmmcluster")
  txdb <- GenomicFeatures::makeTxDbFromGFF(file=gff3, format="gff3") # takes 1 sec, save and load.
  gn <- GenomicFeatures::genes(txdb)
  gn = as.data.frame(gn)
  gn = gn[- which(gn$gene_id %in% ignore_genes),] # remove ignore genes
  # cds and transcripts are both insufficient. preferrably remove certain genes
  # cds <- GenomicFeatures::cds(txdb)
  # cds = GenomicFeatures::transcripts(txdb)
  # cds = GenomicFeatures::transcriptsBy(txdb, by=c("gene", "exon", "cds"), use.names=FALSE)
  # 
  # cds2 = GenomicFeatures::exonicParts(txdb, linked.to.single.gene.only=T)
  # cds2 = as.data.frame(cds2)
  # 
  
  
  
  #annotate data
  dat$genes = ""
  for(r in 1:nrow(dat)){
    if(is.na(regions$ref_start[r]) | is.na(regions$ref_end[r])){next}
    regpos = regions$ref_start[r]:regions$ref_end[r]
    
    genes = ""
    for(g in 1:nrow(gn)){
      gene = gn$gene_id[g]
      gpos = gn$start[g]:gn$end[g]
      condition = regpos %in% gpos
      condition = grep("TRUE", condition)
      condition = length(condition) >= 1
      
      if(condition){
        genes = paste(genes, gene)
      }
      
      # append
      dat$genes[r] = genes
    }
  }
  return(dat)
}
regions = annotate_gene(regions)
write.csv(regions, "regions_annotated.csv",row.names= F)



#----------------------------- 1 genome plot

source("1-circos/circos_function.R") # load function
ignore_genes = c("RL8A", "UL89")
gff3 = system.file("ref/NC_006273.2.gff3", package = "hmmcluster")
ref = system.file("ref/NC_006273.2.fasta", package = "hmmcluster")

plot_circos(gff3,ref,outpng = "1-circos/circos.png",regions = regions, ref_hzyg = "nucdiv-add.csv",ignore_genes = ignore_genes)



#----------------------------- 2 structural and homologous genome + MDS
# will not henikoff weight the large regions as all return 1
dir.create("2-mds")
mds = function(msa){
  pal <- c(
    "Africa" = "#F8766D",
    "Europe" = "#00B81F",
    "Asia"   = "#619CFF"
  )
  dist <- dist.dna( msabin, model = "TN93", pairwise.deletion = T, as.matrix = T )
  
  t.mds <- cmdscale(dist, k = 2)
  t.mds = data.frame(sample = rownames(t.mds),
                     MDS_D1 = t.mds[, 1],
                     MDS_D2 = t.mds[, 2])
  t.mds$continent <-  str_split(t.mds$sample, pattern = "_", simplify = T)[, 1]
  t.mds$country <-   str_split(t.mds$sample, pattern = "_", simplify = T)[, 2]
  t.mds = t.mds[t.mds$continent %in% c("Europe", "Africa", "Asia"),]
  rownames(t.mds) <- NULL
  g = ggplot(t.mds, aes(x = MDS_D1, y = MDS_D2, text = sample, colour = continent, alpha = 0.5, size = 5)) +
    geom_point() +
    theme_classic() +
    guides(alpha=FALSE, size = FALSE)+
    labs(title = paste("HCMV Whole Genome clustering \n by TN93 pairwise distance ")) +
    scale_colour_manual(
      values = pal,
      limits = names(pal)
    ) +
    theme(legend.position = "none")
  
  g1 = g
  ggsave(g1, filename = paste("2-mds/mds-wg.png", sep = ""), device = "png")
  
  
  
  ##--------------- structural
  # generate vector of all structural positions
  struc = c()
  for(reg in regions$region){
    struc = c(struc, regions$start[reg]:regions$end[reg])
  }
  msa.struc = msabin[,struc]
  write.FASTA(msa.struc, "2-mds/msa_struc.fasta")
  dist.struc <- dist.dna( msa.struc, model = "TN93", pairwise.deletion = T, as.matrix = T )
  
  
  t.mds <- cmdscale(dist.struc, k = 2)
  t.mds = data.frame(sample = rownames(t.mds),
                     MDS_D1 = t.mds[, 1],
                     MDS_D2 = t.mds[, 2])
  t.mds$continent <-  str_split(t.mds$sample, pattern = "_", simplify = T)[, 1]
  t.mds$country <-   str_split(t.mds$sample, pattern = "_", simplify = T)[, 2]
  rownames(t.mds) <- NULL
  t.mds = t.mds[t.mds$continent %in% c("Europe", "Africa", "Asia"),]
  g = ggplot(t.mds, aes(x = MDS_D1, y = MDS_D2, text = sample, colour = continent, alpha = 0.5, size = 5)) +
    geom_point() +
    theme_classic() +
    guides(alpha=FALSE, size = FALSE)+
    labs(title = paste("HCMV structural regions clustering \n by TN93 pairwise distance ")) +
    scale_colour_manual(
      values = pal,
      limits = names(pal)
    ) +
    theme(legend.position = "none")
  g2 = g
  ggsave(g, filename = paste("2-mds/mds-struc.png", sep = ""), device = "png")
  
  
  ##--------------- homologous
  # generate vector of all structural positions
  hom = 1:length(msa[1,])
  hom = hom[!(hom %in% struc)]
  
  
  msa.hom = msabin[,hom]
  write.FASTA(msa.hom, "2-mds/msa_hom.fasta")
  dist.hom <- dist.dna( msa.hom, model = "TN93", pairwise.deletion = T, as.matrix = T )
  t.mds <- cmdscale(dist.hom, k = 2)
  t.mds = data.frame(sample = rownames(t.mds),
                     MDS_D1 = t.mds[, 1],
                     MDS_D2 = t.mds[, 2])
  t.mds$continent <-  str_split(t.mds$sample, pattern = "_", simplify = T)[, 1]
  t.mds$country <-   str_split(t.mds$sample, pattern = "_", simplify = T)[, 2]
  rownames(t.mds) <- NULL
  t.mds = t.mds[t.mds$continent %in% c("Europe", "Africa", "Asia"),]
  g = ggplot(t.mds, aes(x = MDS_D1, y = MDS_D2, text = sample, colour = continent, alpha = 0.5, size = 5)) +
    geom_point() +
    theme_classic() +
    guides(alpha=FALSE, size = FALSE)+
    labs(title = paste("HCMV homologous regions clustering \n by TN93 pairwise distance ")) +
    scale_colour_manual(
      values = pal,
      limits = names(pal)
    ) +
    theme(legend.position = "none")
  
  g3 = g
  ggsave(g, filename = paste("2-mds/mds-hom.png", sep = ""), device = "png")
  
  
  #g4 = ggarrange(g1,g2,g3,common.legend = T)
  #ggsave(g4, filename = paste("2-mds/mds-all.png", sep = ""), device = "png",width = 20, height = 20)
  
  leg <- get_legend(g3)
  g4 <- ggarrange(g1, g2, g3,  leg, ncol = 2, nrow = 2)
  ggsave(g4, filename = paste("2-mds/mds-all.png", sep = ""), device = "png",width = 10, height = 10)
  
}
mds(msa)



#----------------------------- 3 phylogenetic tree on all regions, 5 exemplar homologous regions, and homologous genome
dir.create("3-phylo")
system("cp out/3-clean/*.fasta 3-phylo") # copy in region fasta
system("cp hom*.fasta 3-phylo") # copy in example homologous fasta



system(paste0('find 3-phylo/ -name "*.fasta" | parallel iqtree -s  {} -m GTR+G4+FO -bb 1000 -alrt 1000 -asr'),ignore.stdout = T)

# run homologous genome
# command = paste("iqtree -s 2-mds/msa_hom.fasta  -m GTR+G4+FO -bb 1000 -alrt 1000 -asr -nt 12") # allow homologous genome to use up to 12 cores, and do ancestral reconstruction
# system(command, wait = F,ignore.stdout = T)


system("cp 2-mds/*.treefile 3-phylo") # copy in example homologous fasta
source("3.R")



#----------------------------- 4 Recombination
# - do this manually in windows, seems a real hasstle for linux automation

## evidence of recombination in the conserved genome

## evidence of recombination in the structural regions per region


#----------------------------- 5 geo-geno
dir.create("5-geo-geno")

source("5.R")

#----------------------------- 6 geo-fst-scan
dir.create("6-geo-fst")
# gamma st for struc, conserved regions. sample 500 times, at 200bp.

source("6.R")



#----------------------------- 7 Mutual information
dir.create("7-mutinf")
# calculate mutual information for each pair or regions

source("7-mutinf.R")
# todo
# remove the regions that are geographically distributed.


#----------------------------- 8 selection
dir.create("8-selection")




#----------------------------- 10 general stats
# alellicism of positions

ten = list()
ten$alellicism = vector(mode="character", length = ncol(msa))
for(i in 1:ncol(msa)){
  t = as.character(msa[,i])
  t = t[t %in% c("a", "c", "g" ,"t")]
  ten$alellicism[i] = length(unique(t))
}
ten$alellicism = table(ten$alellicism)
print("allellicism")
print(ten$alellicism)


# % of genome that is structural
msa.struc = read.dna("2-mds/msa_struc.fasta", "fasta",as.character = T, as.matrix = T)

print( paste(ncol(msa.struc) / ncol(msa), "% of the genome is structural"))


# size of regions


#----------------------------- 11 extra origins
dir.create("11-origins")
source("11.R")




repeat_regions = c( # likely to be very poor sequenceing using short reads
  TLR = 1:1324, # TLR Terminal Repeat Long
  OAR = 94059:94198, #"OriLyt-associated repea" - two repeats of ~150, good coverage
  IRL = 194344:195667, #IRL Internal Repeat Long
  IRS = 195090:197626, #IRS Internal Repeat Short
  TRS = 233109:235646 # TRS Terminal Repeat Short
)



