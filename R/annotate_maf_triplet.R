#' @title Annotate MAF with triplet context
#'
#' @description Give triple sequence of mutated base with its adjacent bases (-1 and +1) based on the STRAND_VEP column of MAF 
#'
#' @details It gives the reference and alternative alleles and looks for the rows of the data frame based on these values for + strand genes and their 
#' complement alleles rows for - strand genes, then it can look for the adjacent bases in that mutation position. Also, it can look for all the SNVs in 
#' the MAF data frame and provide triple sequences for them.
#'
#' @param maf MAF file (required columns: Reference_Allele, Tumor_Seq_Allele2, STRAND_VEP)
#' @param all_SNVs To give us all the triplet sequences of SNVs and not specifying any specific ref and alt alleles (default is TRUE)
#' @param ref Reference allele
#' @param alt Alternative allele
#' @param projection The genome build projection for the variants you are working with (default is grch37)
#' @param fastaPath Can be a path to a FASTA file
#'
#' @return A data frame with an extra column for triple sequence
#'
#' @import Rsamtools dplyr BSgenome
#' @export
#'
#' @examples
#' annotate_maf_triplet(maf, all_SNVs = FALSE, "C", "T")
#' 

#This function gives triple sequence of provided mutated base
annotate_maf_triplet = function(maf,
                                all_SNVs = TRUE,
                                ref,
                                alt,
                                projection = "grch37",
                                fastaPath,
                                bsgenome_name,
                                pyrimidine_collapse = FALSE){
  genome = ""
  bsgenome_loaded = FALSE
  if (projection == "grch37") {
    maf$Chromosome <- gsub("chr", "", maf$Chromosome)
  } else {
    # If there is a mix of prefixed and non-prefixed options
    maf$Chromosome <- gsub("chr", "", maf$Chromosome) 
    maf$Chromosome <- paste0("chr", maf$Chromosome)
  }
  # If there is no fastaPath, it will read it from config key 
  # Based on the projection the fasta file which will be loaded is different
  if (missing(fastaPath)){
    
    
    base <- check_config_value(config::get("repo_base"))
    fastaPath <- paste0(
      base,
      "ref/lcr-modules-references-STABLE/genomes/",
      projection,
      "/genome_fasta/genome.fa"
    )
    if(!file.exists(fastaPath)){
      #try BSgenome
      installed = installed.genomes()
      if(!missing(bsgenome_name)){
        if(bsgenome_name %in% installed){
          genome = getBSgenome(bsgenome_name)
          
          bsgenome_loaded = TRUE
        }
      }else{
        print(installed)
        stop("Local Fasta file cannot be found. Supply a path to it with fastaPath or use one of the installed BSGenome options using the bsgenome_name parameter")
      }
    }
  }
  # It checks for the presence of a local fastaPath
  if(!bsgenome_loaded){
    # Create a reference to an indexed fasta file.
    if (!file.exists(fastaPath)) {
      stop("Failed to find the fasta file")
    }
    fasta = Rsamtools::FaFile(file = fastaPath)
  }

  # Store the complement of ref and alt alleles
  complement <- c(
    'A'= 'T',
    'T'= 'A',
    'C'= 'G',
    'G'= 'C'
  )  
  if(pyrimidine_collapse){
    maf = mutate(maf,mutation_strand=ifelse(Reference_Allele %in% c("A","G"),"-","+")) %>%
      mutate(mutation = ifelse(Reference_Allele %in% c("A","G"),
                               paste0(complement[Reference_Allele],">",complement[Tumor_Seq_Allele2]),
             paste0(Reference_Allele,">",Tumor_Seq_Allele2)))
  }else{
    maf = mutate(maf,mutation_strand="+")
  }
  
  CompRef = complement[ref]
  CompAlt = complement[alt]
  # Provide triple sequence for all the SNVs
  if (all_SNVs == TRUE){
    if(bsgenome_loaded){
      sequences <- maf %>%
        dplyr::mutate(
          seq = ifelse(
            (nchar(maf$Reference_Allele) == 1 &
               nchar(maf$Tumor_Seq_Allele2) == 1
            ),
            as.character(getSeq(
              genome,
              maf$Chromosome,
              start = maf$Start_Position-1,
              end = maf$Start_Position + 1,
              strand = mutation_strand
              )),
            "NA"
          )
        )
      
    }else{

      #TODO: if we want this to work with a local fasta we need to ensure the - strand is returned. I couldn't test it on my Laptop
      sequences <- maf %>%
        dplyr::mutate(
          seq = ifelse(
            (nchar(maf$Reference_Allele) == 1 &
               nchar(maf$Tumor_Seq_Allele2) == 1
            ),
            Rsamtools::getSeq(
              fasta,
              GenomicRanges::GRanges(
                maf$Chromosome,
                IRanges(
                  start = maf$Start_Position - 1,
                  end = maf$Start_Position + 1
                )
              )
            ),
            "NA"
          )
        )
      
    }
    
  }else{
    # Provide triple sequence of + strand and reverse complement of - strand   
    # Mutations on + strand with chosen ref and alt alleles
    # Mutations on - strand with complement ref and alt alleles
    if (!"STRAND_VEP" %in% colnames(maf)){
      stop("STRAND_VEP column is not in the MAF")
    }
    sequences <- maf %>%
      dplyr::mutate(
        seq = ifelse(
          (maf$STRAND_VEP == "1" &
             maf$Reference_Allele == ref &
             maf$Tumor_Seq_Allele2 == alt
          )|(
            maf$STRAND_VEP == "-1" &
              maf$Reference_Allele == CompRef &
              maf$Tumor_Seq_Allele2 == CompAlt
          ),
          Rsamtools::getSeq(
            fasta,
            GenomicRanges::GRanges(
              maf$Chromosome,
              IRanges(
                start = maf$Start_Position - 1,
                end = maf$Start_Position + 1
              )
            )
          ),
          "NA"
        )
      )  
    
  }
  return(sequences)
}
