#' Handles input file
#'
#' This function handles a vcf, or varscan tab variant file, or pases fasta files
#' onto handle_fasta().
#' Returns an intermediate data.frame, which contains variant information to be
#' annotated by annotate_variants().
#'
#' @param f.infile Path to the input file
#' @param global Package object for consistent runtime variables
#' @return An intermediate data.frame
read_input <- function(f.infile, global){
  #takes .vcf .tab .fasta inputs and returns a standard dataframe for resistance calling
  f.infile <- as.character(f.infile)
  ### tab ###
  if(tools::file_ext(f.infile) == "tab"){
    #if file appears as a varscan tab file
    # delimit
    tab.dat <- read.table(file = f.infile, header = T, as.is = T, sep = "\t")
    out <- read_varscan_data(tab.dat)
  }
  ### vcf ###
  else if(tools::file_ext(f.infile) == "vcf"){
    
    text <- readLines(f.infile)
    start <- grep('chrom',ignore.case = T, text)
    vcf = read.delim(f.infile, sep = "\t", as.is = T, skip = start - 1)
    
    # steves vcf files are occasionally ...dodge
    vcf = vcf[grepl(pattern = ".{15,100}",vcf[,10]),] # where the final info column is actually legitimate
    vcf = vcf[!grepl(pattern = "\\.",vcf[,5]),] # remove any positions with alt variant of . - this seems obvious. remvoes on variants
    
    
    
    # get format column
    t.col = grep(pattern = "FORMATa",x = names(vcf))
    
    if(length(vcf[,9]) > 0){ # if has a format column & genotype column, split to extract ref.count, var.count per position
      vcf.num_format = as.numeric(length(unlist(strsplit(vcf[1,9], split=":"))))
      t.1 <- as.data.frame(matrix(unlist(strsplit(vcf[,10], split=":")), ncol=vcf.num_format, byrow="T"), stringsAsFactors=F)
      colnames(t.1) <- unlist(strsplit(vcf[1,9], split=":"))
      
    }else{
      stop("Check your variant call file has genotypic information, it mght be there but isnt in the standard format!")
    }
    
    for(i in 1:nrow(vcf)){#clean up vcf indel format to be as in varscan tab
      ref = vcf$REF[i]
      var = vcf$ALT[i]
      if(nchar(ref) > 1){#if deletion
        out.ref = var
        out.var = ref
        substr(out.var, 1, 1) <- "-"
        vcf$REF[i] = out.ref
        vcf$ALT[i] = out.var
      }
      if(nchar(var) > 1){#if insertion
        out.ref = ref
        out.var = var
        substr(out.var, 1, 1) <- "+"
        vcf$REF[i] = out.ref
        vcf$ALT[i] = out.var
      }
    }
    
    t.vcf <- data.frame(Position = vcf$POS,
                        Ref = vcf$REF,
                        Var = vcf$ALT,
                        Ref.count = t.1$RD,
                        Var.count = t.1$AD,
                        VarFreq = t.1$FREQ,
                        Sample = "single run",
                        stringsAsFactors = F)
    out <- t.vcf
    
    
  }
  ### fasta ###
  else if(tools::file_ext(f.infile) %in% c("fa", "fasta", "fas")){
    # in each case output is a vcf file, which then gets processed as above into the out data structure.
    
    #writes vcf to a known location accessible by global$dir variable
    vcf_file = handle_fasta(fasta_in = f.infile, fasta_out = paste(global$dir, "out.fasta", sep = ""), fasta_ref = global$path_fasta_file) 
    
    text <- readLines(vcf_file)
    start <- grep('chrom',ignore.case = T, text)
    vcf = read.delim(vcf_file, sep = "\t", as.is = T, skip = start - 1)
    
    # if has a format column & genotype column, split to extract ref.count, var.count per position
    for(i in 1:nrow(vcf)){#clean up vcf indel format to be as in varscan tab
      ref = vcf$REF[i]
      var = vcf$ALT[i]
      if(nchar(ref) > 1){#if deletion
        out.ref = var
        out.var = ref
        substr(out.var, 1, 1) <- "-"
        vcf$REF[i] = out.ref
        vcf$ALT[i] = out.var
      }
      if(nchar(var) > 1){#if insertion
        out.ref = ref
        out.var = var
        substr(out.var, 1, 1) <- "+"
        vcf$REF[i] = out.ref
        vcf$ALT[i] = out.var
      }
    }
    
    t.vcf <- data.frame(Position = vcf$POS,
                        Ref = vcf$REF,
                        Var = vcf$ALT,
                        Ref.count = vcf[,10], #diff from vcf proc
                        Var.count = vcf[,11], # diff from vcf proc
                        VarFreq = "100%",
                        Sample = "single run",
                        stringsAsFactors = F)
    out <- t.vcf
    
    
    
    
  }else{
    stop("Check your variant call file is in .tab, .vcf format \n or check your fasta file has the .fa, .fas or .fasta extension")
    
  }
  
  return(out)
  #remember to update read_Varscan input functions to read a dataframe not a file location
}
