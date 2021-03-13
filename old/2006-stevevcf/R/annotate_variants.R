#' Annotates variants
#'
#' Adds genome annotation to the intermediate resistance data.frame
#'
#' @param f.dat intermediate data.frame
#' @param global Package object for consistent runtime variables
#' @return intermediate data.frame with genome level annotation
#' 

annotate_variants <- function(f.dat,global){
  toannotate <- f.dat
  check <- IRanges(start=toannotate$Position, end=toannotate$Position, width=1)
  gr <- GRanges(seqnames = global$genome ,ranges=check)
  values(gr) <- DataFrame(id = toannotate$Sample, freq = toannotate$VarFreq, RefCount= toannotate$Ref.count, VarCount= toannotate$Var.count, VarAllele=toannotate$Var)
  varallele <- DNAStringSet(toannotate$Var)
  txdb <- makeTxDbFromGFF(file=global$path_gff3_file, format="gff3")
  #select(txdb, columns=columns(txdb), keys=keys(txdb), keytype=c("GENEID"))
  gn <- genes(txdb)
  cds(txdb)
  ##variant data
  vcf <- gr
  #seqlevels need to match with txdb
  #vcf <- gr
  seqlevels(txdb) <- global$genome
  seqlengths(vcf) <- seqlengths(txdb)[names(seqlengths(vcf))]
  isCircular(vcf) <- isCircular(txdb)[names(seqlengths(vcf))]
  #genome sequence from fasta file
  fa <- FaFile(global$path_fasta_file)
  codvar <- suppressWarnings(locateVariants(vcf, txdb, CodingVariants()))
  coding <- suppressWarnings(predictCoding(vcf, varAllele = varallele,txdb, seqSource=fa , ignore.strand=FALSE))
  txids <- values(transcripts(txdb))$tx_name
  names(txids) <- values(transcripts(txdb))$tx_id
  coding_df <- as.data.frame(coding)
  #coding_df <- DataFrame(coding)
  coding_df$aachange <- paste(coding_df$REFAA,coding_df$PROTEINLOC,coding_df$VARAA,sep="")
  coding_df$change <- paste(coding_df$GENEID,coding_df$aachange,sep="_")
  return(coding_df)
}
