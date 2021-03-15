# for each patient in the table d2 - only keep mutations that occur greater than 2* at5% or more. -> oper patient

infile = "data/2103-SuppTable_RawVariants.csv"
outdir = "analysis/2103-perpatientvarmuts/"

library(stringr)


df = read.csv(infile)
df[,c(1,3,4,5)] <- NULL

colnames(df) = str_remove_all(colnames(df), "X")


cols = str_split(colnames(df), "_",simplify = T)[,1]

patients = cols[2:75]

patients = unique(patients)

for(patient in patients){

  #patient = patients[1]
  
  wcol = grep(patient, cols) # which cols

  
  tdf = df[,c(1,wcol)]
  
  
  
  # filter as per steves wishes
  # keep a mutations if it is above 5% in at lest 2 tp's.
  keep = c()
  for(r in 1:nrow(tdf)){
    t = as.numeric(tdf[r,-1])
    t = length(which(t > 5))
    if(t >= 2){ keep = c(keep, r)} # vector or rows to keep
  }
  out = tdf[keep,]
  
  write.csv(out, paste0("analysis/2103-perpatientvarmuts/", patient, ".csv"))
  
  
  
}