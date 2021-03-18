# calculate dnds for gag,pol,env
library(seqinr)
library(Biostrings)
library(dplyr)
outdir = "analysis/7_dnds/"
msa_file = "old/2102-mds/fig1/consensus_filt.fasta"

msa = ape::read.dna(msa_file,format = "fasta", as.matrix = T, as.character = T)

regions = c("gag" = (792-3):2279, "pol" = 2078:5083, "env" = 6230:8812)

#-------------------- run manually
gene = "env"
t = msa[,6230:8812] # include Methionine alter this line
ape::write.dna(t,paste0(outdir, gene,"-dna.fasta"), "fasta")


# manually edit the alignment to be by translation - use aliview,better than anything automated :(

t = ape::read.dna(paste0(outdir, gene,"-dna-man.fasta"),format = "fasta",as.character = T, as.matrix = T)
starts = seq(1,ncol(t),3)
df = data.frame(starts = starts)
df$kaks = 0
for(i in 1:length(starts)){
  start =  starts[i]
  t2 = t[,start:(start+2)]
  ape::write.dna(t2, "temp.fasta",format = "fasta")
  t3 = seqinr::read.alignment("temp.fasta",format = "fasta")
  t4 = seqinr::kaks(t3)
  df$kaks[i] = mean(t4$ka) / mean(t4$ks)
  print(mean(t4$ka) / mean(t4$ks))  
}
file.remove("temp.fasta")
df = df[df$kaks != 0,]
df = df[df$kaks != 1,]
write.csv(df, paste0(outdir, gene,".csv"))
plot(df)
png(paste0(outdir,gene,".png"))
