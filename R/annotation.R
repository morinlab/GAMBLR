#functions for custom annotation of mutations such as SV, CNV

#some data frames specifying the locations of oncogenes and common partners are included with this package


#' Retrieve data or use supplied data frame of genes to get/annotate MAF data indicating which rows are putative driver mutations
#'
#' @param maf_df Optional data frame of MAF-format mutations (default is to retrieve automatically)
#' @param lymphoma_type Optional keyword to find genes for annotation (e.g. "BL","DLBCL","FL") is not specifying driver genes
#' @param driver_genes Optional vector of Hugo_Symbol of genes for coding mutations
#' @param include_noncoding Optional named vector of genes listing the Variant_Classification of noncoding mutations to retain
#' @param noncoding_regions Optional named vector of regions to use to further restrict noncoding mutations per gene
#'
#' @return
#' @export
#'
#' @examples
#' driver_ssm = annotate_driver_ssm(include_noncoding=c("NFKBIZ"="3'UTR"),noncoding_regions=c("NFKBIZ"="chr3:101578206-101578365"))
annotate_driver_ssm = function(maf_df,lymphoma_type,driver_genes,
                               include_noncoding=c("NFKBIZ"="3'UTR","HNRNPH1"="Intron"),
                               noncoding_regions=c("NFKBIZ"="chr3:101578206-101578365","HNRNPH1"="chr5:179,045,946-179,046,683")){
  coding_vc = c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del","In_Frame_Ins","Missense_Mutation","Nonsense_Mutation","Nonstop_Mutation","Splice_Region","Splice_Site","Targeted_Region","Translation_Start_Site")
  #get the gene list if missing
  if(missing(driver_genes)){
    driver_genes = lymphoma_genes[which(lymphoma_genes[[lymphoma_type]]==TRUE),] %>% pull(Gene)
    ngene = length(driver_genes)
    message(paste("using", ngene, "genes"))

  }
  if(missing(maf_df)){
    #use the database to get the coding and noncoding variants for the user
    kept_ssm = get_coding_ssm() %>% dplyr::filter(Variant_Classification %in% coding_vc) %>%
      dplyr::filter(Hugo_Symbol %in% driver_genes)
    message("starting with:")
    print(dim(kept_ssm))
    #this doesn't yet have any of our noncoding variants
    for(gene in names(include_noncoding)){
      message(paste("adding",unname(include_noncoding[gene]), "for", gene))
      nc_ssm = get_ssm_by_gene(gene_symbol = "NFKBIZ") %>%
        dplyr::filter(Variant_Classification == unname(include_noncoding[gene]))
      if(!is.na(noncoding_regions[gene])){
        #also restrict to coordinates
        chunks = region_to_chunks(noncoding_regions[gene])
        nc_ssm = dplyr::filter(nc_ssm,Start_Position >= chunks$start & Start_Position <= chunks$end)
        message(gene)
        print(dim(nc_ssm))
      }
      kept_ssm = bind_rows(kept_ssm,nc_ssm)
      print(dim(kept_ssm))
    }
  }else{
    #this all assumes the user has left the noncoding variants in their MAF
    kept_ssm = maf_df %>% dplyr::filter(Variant_Classification %in% coding_vc) %>%
      dplyr::filter(Hugo_Symbol %in% driver_genes)
    for(gene in names(include_noncoding)){
      message(paste("adding",unname(include_noncoding[gene]), "for", gene))
      nc_ssm = maf_df %>% dplyr::filter(Variant_Classification == unname(include_noncoding[gene]))
      if(!is.na(noncoding_regions[gene])){
        #also restrict to coordinates
        chunks = region_to_chunks(noncoding_regions[gene])
        nc_ssm = dplyr::filter(nc_ssm,Start_Position >= chunks$start & Start_Position <= chunks$end)
      }
      kept_ssm = bind_rows(kept_ssm,nc_ssm)
      print(dim(kept_ssm))
    }

  }
  return(kept_ssm)


}


#' Annotate a data frame of SV breakpoints after retrieval from the database
#'
#' @param sv_data A data frame of SVs. This should be the output of get_manta_sv. If you aren't using the database backend you can supply your own data frame in the format show below.
#' Most of this data is directly from the bedpe files that are obtained by converting the Manta outputs from VCF.
#' Only the columns for both chromosomes, coordinates and strand plus SOMATIC_SCORE and tumour_sample_id are absolutely required
#'  CHROM_A  START_A    END_A CHROM_B  START_B    END_B NAME SOMATIC_SCORE STRAND_A STRAND_B TYPE FILTER VAF_tumour VAF_normal DP_tumour DP_normal tumour_sample_id normal_sample_id pair_status
#'   1  1556541  1556547       1  1556664  1556670    .            40        -        -  BND   PASS      0.145          0        55        73  00-14595_tumorA  00-14595_normal     matched
#' @param partner_bed Optional bed-format data frame to use for annotating oncogene partners (e.g. enhancers). required columns are: chrom,start,end,gene
#' @param with_chr_prefix Optionally request that chromosome names are returned with a chr prefix
#' @param collapse_redundant Remove reciprocal events and only return one per event
#'
#' @return A data frame with annotated SVs (gene symbol and entrez ID)
#' @export
#' @import tidyverse
#' @import data.table
#'
#' @examples
#' # Basic usage
#' sv_df = get_manta_sv()
#' annotated_sv = annotate_sv(sv_df)
annotate_sv = function(sv_data,partner_bed,with_chr_prefix=FALSE,collapse_redundant=FALSE,return_as="bedpe"){
  bedpe1 = sv_data %>% dplyr::select("CHROM_A","START_A","END_A","tumour_sample_id","SOMATIC_SCORE","STRAND_A")
  bedpe2 = sv_data %>% dplyr::select("CHROM_B","START_B","END_B","tumour_sample_id","SOMATIC_SCORE","STRAND_B")

  colnames(bedpe1)= c("chrom","start","end","tumour_sample_id","score","strand1")
  colnames(bedpe2)= c("chrom","start","end","tumour_sample_id","score","strand2")
  suppressWarnings({
    if(grepl("chr",bedpe1$chrom)){
      bedpe1 = dplyr::mutate(bedpe1,chrom = str_replace(chrom,"chr",""))
      bedpe2 = dplyr::mutate(bedpe2,chrom = str_replace(chrom,"chr",""))
    }
  })
  if(missing(partner_bed)){
    ig_regions = grch37_partners
  }else{
    ig_regions = partner_bed
    if(!"entrez" %in% colnames(ig_regions)){
      ig_regions$entrez = 0
    }
  }
  oncogene_regions = grch37_oncogene
  y = data.table::as.data.table(oncogene_regions)

  data.table::setkey(y, chrom, start, end)
  #use foverlaps to get oncogene SVs

  a = data.table::as.data.table(bedpe1)
  a.onco = data.table::foverlaps(a, y, type="any",mult="first") #oncogene-annotated bedpe for the first breakpoints
  b = data.table::as.data.table(bedpe2)
  b.onco = data.table::foverlaps(b, y, type="any",mult="first") #oncogene-annotated bedpe for the first breakpoints

  #insist oncogene breakpoints are anchored in an IG or superenhancer region (currently just IG or BCL6)
  #other end of breakpoint

  a.onco.break = a.onco[which(!is.na(a.onco$start)),c("chrom","i.start","i.end","tumour_sample_id","gene","entrez","score","strand1")]
  b.onco.break = b.onco[which(!is.na(b.onco$start)),c("chrom","i.start","i.end","tumour_sample_id","gene","entrez","score","strand2")]

  a.partner = b[which(!is.na(a.onco$start)),]
  b.partner = a[which(!is.na(b.onco$start)),]

  y = data.table::as.data.table(ig_regions)

  data.table::setkey(y, chrom, start, end)

  a.ig = data.table::foverlaps(a.partner, y, type="any",mult="first")

  b.ig = data.table::foverlaps(b.partner, y, type="any",mult="first")

  a.ig = a.ig[,c("chrom","i.start","i.end","strand2","gene")]
  b.ig = b.ig[,c("chrom","i.start","i.end","strand1","gene")]
  a.annotated.both = cbind(a.onco.break,a.ig)
  colnames(a.annotated.both) = c("chrom1","start1","end1","tumour_sample_id","gene","entrez","score","strand1","chrom2","start2","end2","strand2","partner")

  b.annotated.both = cbind(b.onco.break,b.ig)
  colnames(b.annotated.both) = c("chrom2","start2","end2","tumour_sample_id","gene","entrez","score","strand2","chrom1","start1","end1","strand1","partner")

  all.annotated = rbind(a.annotated.both,b.annotated.both)
  all.annotated$fusion = dplyr::pull(tidyr::unite(all.annotated,fusion,partner,gene,sep="-"),fusion)
  all.annotated = dplyr::filter(all.annotated,!fusion %in% c("BCL6-BCL6","CIITA-CIITA","FOXP1-FOXP1"))

  #TODO: need a better system for cataloguing and using these but this works for our current data (hg19 coordinates)
  blacklist = c(60565248,30303126,187728894,101357565,101359747,161734970,69400840,65217851,187728889,188305164)

  all.annotated  = dplyr::filter(all.annotated,!start1 %in% blacklist)
  all.annotated  = dplyr::filter(all.annotated,!start2 %in% blacklist)

  if(return_as=="bedpe"){
    all.annotated$name= "."
    all.annotated = dplyr::select(all.annotated,chrom1,start1,end1,chrom2,start2,end2,name,score,strand1,strand2,tumour_sample_id,gene,partner,fusion)
  }else if(return_as == "bed"){
    #lose the linkage but add a name that somewhat retains it
    if(!grepl("chr",all.annotated$chrom1)){
      all.annotated = all.annotated %>% dplyr::mutate(chrom1 = paste0("chr",chrom1))
      all.annotated = all.annotated %>% dplyr::mutate(chrom2 = paste0("chr",chrom2))
    }
    bed1 = dplyr::mutate(all.annotated,name=paste(tumour_sample_id,fusion,sep="_")) %>% dplyr::select(chrom1,start1,end1,name,score,strand1)
    bed2 = dplyr::mutate(all.annotated,name=paste(tumour_sample_id,fusion,sep="_")) %>% dplyr::select(chrom2,start2,end2,name,score,strand2)
    colnames(bed1)=c("chrom","start","end","name","score","strand")
    colnames(bed2)=c("chrom","start","end","name","score","strand")
    return(dplyr::arrange(rbind(bed1,bed2),name))
  }else{
        if(collapse_redundant){
      all.annotated = dplyr::distinct(all.annotated,tumour_sample_id,fusion,.keep_all = TRUE)
    }

  }
  if(with_chr_prefix){
    #add the prefix if necessary
    if(!grepl("chr",all.annotated$chrom1)){
      all.annotated = all.annotated %>% dplyr::mutate(chrom1 = paste0("chr",chrom1))
      all.annotated = all.annotated %>% dplyr::mutate(chrom2 = paste0("chr",chrom2))
    }
  }
  return(all.annotated)
}

