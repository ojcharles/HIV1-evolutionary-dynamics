# March 2021 - this generates the d2 table, where each column is a patient, by mutation and value is frequency
# reads a csv per (patient, timepoint), and returns a table with columns

library(stringr)
library(reshape2)
library(readr)


files = list.files("data/vcf-mar21-trust/",full.names = T) # vcf-csv1 if the respository that steve agreed is the correct vcf files
files = files[grepl(".csv", files)]


# for each vcf.csv append the data to a large composite dataframe.
out = data.frame()
for(i in 1:length(files)){

  #i = 1
  file = files[i]
  
  d = read.csv(file,stringsAsFactors = F)
  d = d[,c("GENEID", "aachange", "freq", "CONSEQUENCE")]
  pat_tp = str_split(basename(file), ".vcf",simplify = T)[1]
  d2 = cbind(patient_tp = pat_tp, d)
  
  out = rbind(out, d2)

}

#remove indels
keep = grepl("[A-Z]{1}[0-9]{1,4}[A-Z]{1}",out$aachange)
out = out[keep,]
out$freq = parse_number(out$freq)


### remove mutations that are consistently under 5% or consistenelty over 95%
out$pat = str_split(out$patient_tp, "_", simplify = T)[,1]
out$tp = str_split(out$patient_tp, "_", simplify = T)[,2]
out$protmut = paste0(out$GENEID, "_", out$aachange)

out2 = data.frame()
for(patient in unique(out$pat)){
  t = out[out$pat == patient,]
  keep = c("")
  #mut = "env_A532A"
  for(mut in unique(t$protmut)){
    t1 = t[t$protmut == mut,]
    if(nrow(t1) == 1){keep = c(keep, mut)} # sequence failure, keep
    #else if(min(t1$freq) > 95){next}
    #else if(max(t1$freq) < 5){next}
    else{keep = c(keep,mut)
    }
  }
  t = t[t$protmut %in% keep,]
  if(nrow(out2) ==0){
    out2 = t
  }else{
    out2 = rbind(out2, t)
  }
}


#dcast
df = dcast(out2, GENEID + aachange + CONSEQUENCE~ patient_tp, value.var = "freq",fun.aggregate = max)
df[df == -Inf] = 0
df$change = paste0(df$GENEID, "_",df$aachange)

# res merge
res = read.table("data/resmuts2.tsv",sep = "\t", stringsAsFactors = F)
res$change = paste0(res$V1, "_", res$V2)
res$phenotype = "Resistant"

df = merge(df, res, by = "change", all.x = T )
phenotype = df$phenotype
df$phenotype <- NULL
df = cbind(phenotype, df)
df$V1 <- NULL
df$V2 <- NULL
write.table(df, "analysis/2_patientvcf2data/d2a.tsv",row.names = F,sep = "\t")







