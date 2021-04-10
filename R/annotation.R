require("dbplyr")
require("tidyverse")
require("data.table")

#functions for custom annotation of mutations such as SV, CNV

#some data frames specifying the locations of oncogenes and common partners are included with this package

#' Annotate a data frame of SV breakpoints after retrieval from the database
#'
#' @param sv_data A data frame of SVs. This should be the output of get_manta_sv.
#' @param genome_build Which reference genome these SVs are from
#'
#' @return A data frame with annotated SVs (gene symbol and entrez ID)
#' @export 
#'
#' @examples 
#' # Basic usage
#' annotated_sv = annotate_sv(sv_df)
annotate_sv = function(sv_data,genome_build="grch37"){
  bedpe1 = sv_data %>% select("CHROM_A","START_A","END_A","tumour_sample_id")
  bedpe2 = sv_data %>% select("CHROM_B","START_B","END_B","tumour_sample_id")
  
  colnames(bedpe1)[c(1:3)]= c("chrom","start","end")
  colnames(bedpe2)[c(1:3)]= c("chrom","start","end")
  if(genome_build == "grch37" || genome_build == "hg19"){
    oncogene_regions = grch37_oncogene
    ig_regions = grch37_partners
  }
  y = as.data.table(oncogene_regions)
  
  setkey(y, chrom, start, end)
  #use foverlaps to get oncogene SVs
  
  a = as.data.table(bedpe1)
  a.onco = foverlaps(a, y, type="any") #oncogene-annotated bedpe for the first breakpoints
  b = as.data.table(bedpe2)
  b.onco = foverlaps(b, y, type="any") #oncogene-annotated bedpe for the first breakpoints
  
  #insist oncogene breakpoints are anchored in an IG or superenhancer region (currently just IG or BCL6)
  #other end of breakpoint
  
  a.onco.break = a.onco[which(!is.na(a.onco$start)),c("chrom","i.start","i.end","tumour_sample_id","gene","entrez")]
  b.onco.break = b.onco[which(!is.na(b.onco$start)),c("chrom","i.start","i.end","tumour_sample_id","gene","entrez")]
  
  a.partner = b[which(!is.na(a.onco$start)),]
  b.partner = a[which(!is.na(b.onco$start)),]
  
  y = as.data.table(ig_regions)
  #y = as.data.table(ig_regions.hg19)
  setkey(y, chrom, start, end)
  
  a.ig = foverlaps(a.partner, y, type="any",mult="first") 
  
  b.ig = foverlaps(b.partner, y, type="any",mult="first") 
  
  a.ig = a.ig[,c("chrom","i.start","i.end","gene")]
  b.ig = b.ig[,c("chrom","i.start","i.end","gene")]

  a.annotated.both = cbind(a.onco.break,a.ig)
  colnames(a.annotated.both) = c("chrom1","start1","end1","tumour_sample_id","gene","entrez","chrom2","start2","end2","partner")
  b.annotated.both = cbind(b.onco.break,b.ig)
  colnames(b.annotated.both) = c("chrom1","start1","end1","tumour_sample_id","gene","entrez","chrom2","start2","end2","partner")
  
  all.annotated = rbind(a.annotated.both,b.annotated.both)
  all.annotated$fusion = pull(unite(all.annotated,fusion,partner,gene,sep="-"),fusion)
  all.annotated = dplyr::filter(all.annotated,fusion != "BCL6-BCL6")
  
  #TODO: need a better system for cataloguing and using these
  blacklist = c(60565248,30303126,187728894,101357565,101359747,161734970,69400840,65217851,187728889)
  
  all.annotated  = dplyr::filter(all.annotated,!start1 %in% blacklist)
  all.annotated  = dplyr::filter(all.annotated,!start2 %in% blacklist)
  all.annotated = distinct(all.annotated,tumour_sample_id,fusion,.keep_all = TRUE)
  return(all.annotated)
}
