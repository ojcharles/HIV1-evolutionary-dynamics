# calculate dnds for gag,pol,env
library(seqinr)
library(Biostrings)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(stringr)

### variables
outdir = "analysis/7_dnds/"
msa_file = "old/2102-mds/fig1/consensus_filt.fasta"
region = data.frame(gene = c("gag","protease","reverse transcriptase","RNase","integrase","pol","env"),
           start = c((792-3),2237,2534,3854,4214,2078,6230),
           end = c(2279,2533,3853,4213,5080,5083,8812))
run_from_scratch = F
### end



msa = ape::read.dna(msa_file,format = "fasta", as.matrix = T, as.character = T)
#-------------------- gene -> dnds and plot. manually check each alignment!
for(r in 1:nrow(region)){ # for each row
  gene = region$gene[r]
  if(run_from_scratch == T){
    t = msa[,(region$start[r]):(region$end[r])] # include Methionine alter this line
    ape::write.dna(t,paste0(outdir, gene,"-dna.fasta"), "fasta") # write msa to manual editing
    
  
    t = ape::read.dna(paste0(outdir, gene,"-dna-man.fasta"),format = "fasta",as.character = T, as.matrix = T) # manually edit and use this
    starts = seq(1,ncol(t),3)
    df = data.frame(starts = starts)
    df$kaks = 0
    for(i in 1:length(starts)){ # calculate the dnds per DNA position
      start =  starts[i]
      t2 = t[,start:(start+2)]
      ape::write.dna(t2, "temp.fasta",format = "fasta")
      t3 = seqinr::read.alignment("temp.fasta",format = "fasta")
      t4 = seqinr::kaks(t3)
      df$kaks[i] = mean(t4$ka) / mean(t4$ks)
      print(mean(t4$ka) / mean(t4$ks))  
    }
    file.remove("temp.fasta")
    df = df[df$kaks != 0,] # remove 0 and 1, uninteresting
    df = df[df$kaks != 1,]
    write.csv(df, paste0(outdir, gene,".csv"),row.names = F)
  }
  
  df = read.csv(paste0(outdir, gene,".csv"))
  df = df[!is.na(df$starts),] # remove random na
  ### read in resistance mutations and merge
  drm = read.delim("data/resmuts2.tsv",sep = "\t")
  drm$aaloc = str_extract(drm$mutation,"[0-9]{1,4}")
  drm$gene[drm$gene == "rt"] <- "reverse transcriptase"
  drm = drm[drm$gene == gene,]
  # 
  # DRms concatenate all mutation by aaloc
  drm = drm %>%
    group_by(aaloc) %>%
    mutate(mut_by_loc = paste0(mutation, collapse = " "))
  
  
  # identify mutations that are resistant and assign a name
  df$aaloc = as.numeric(df$starts) / 3 ;  df$aaloc = ceiling(df$aaloc) # get aaloc from dnapos
  df$drm = "no"
  for(i in 1:nrow(df)){ # for each row if theres an entry in drm mark is as resistant
    aaloc = df[i,]
    aaloc = merge(aaloc,drm, by = "aaloc")
    if(nrow(aaloc) > 0){
      df$drm[i] = aaloc$mutation
    }
  }
  
  # alter mutations to be gene_loc*
  df$drm = str_extract(df$drm, "[A-Z]{1}[0-9]{1,4}")
  df$drm[grepl("[0-9]",df$drm)] <- paste0(df$drm[grepl("[0-9]",df$drm)],"*")
  
  
  if(gene == "gag"){
    # add pi associated
    pi = read.csv("data/gag-pi-associated-sites.csv")[1:53,1:2]
    pi$aaloc = pi$position
    # identify mutations that are resistant and assign a name
    for(i in 1:nrow(df)){ # for each row if theres an entry in drm mark is as resistant
      aaloc = df[i,]
      aaloc = merge(aaloc,pi, by = "aaloc")
      if(nrow(aaloc) > 0){
        df$drm[i] = aaloc[1,5]
      }
    }
  }
  

  ggplot(df)+
    geom_point(aes(x = starts, y = kaks))+
    geom_point(data = df[df$drm != "no",], mapping = aes(x = starts, y = kaks),colour = "red") +
    ggrepel::geom_text_repel(data = df[df$drm != "no",], mapping = aes(x = starts, y = kaks, label = drm),colour = "red", size = 3) +
    theme_classic() +
    scale_x_continuous (breaks = c(1,250,500,750,1000,1250,1500,1750,2000,2250,2500,3000,3250,3500)) +
    geom_hline(aes(yintercept = 1), colour = "red", linetype = "dashed") +
    labs(x=paste(gene,"nucleotide position"), y="dN/dS")
  ggsave(paste0(outdir,gene,".png"))
}