#adding coding_vc to global enviroment
coding_vc = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Region", "Splice_Site", "Targeted_Region", "Translation_Start_Site")


#' @title Annotate SSM with Blacklists
#'
#' @description Annotate and auto-drop a MAF data frame with existing blacklists.
#'
#' @details Annotate and auto-drop a MAF data frame with existing blacklists to remove variants that would be dropped during the merge process.
#' This function returns a MAF format data frame with two new columns, indicating the number of occurrences of each variant in the two blacklists.
#' Note that there are a collection of parameters to this function to improve flexibility for many applications,
#' such as `return_blacklist` (returns the used blacklist to the vector given the function, or printed to the terminal if blank).
#' For returning variants that would be dropped, one can specify `invert = TRUE`, please use with caution, this is most likely the opposite of what you want from this function.
#' Lastly, the minimum count from one of the blacklists to drop a variant is specified with `drop_threshold = 4`.
#' This function also conveniently lets you know how many variants that were dropped in the annotation process.
#'
#' @param mutations_df A data frame with mutation data.
#' @param seq_type The seq_type of your mutations if you prefer to apply only the corresponding blacklist. More than one seq_type can be specified as a vector if desired. This parameter is required.
#' @param tool_name The tool or pipeline that generated the files (should be the same for all).
#' @param tool_version The version of the tool specified under `tool_name`.
#' @param annotator_name Name of annotator, default is "vcf2maf".
#' @param annotator_version Version of annotator specified under `annotator_name`.
#' @param genome_build The genome build projection for the variants you are working with (default is grch37).
#' @param project_base Optional: A full path to the directory that your blacklist_file_pattern is relative to.
#' @param blacklist_file_template Optional: A string that contains the relative path to your blacklist file from after the project_base (i.e. results) with any wildcards surrounded with curly braces.
#' @param drop_threshold The minimum count from one of the blacklists to drop a variant.
#' @param return_blacklist Boolean parameter for returning the blacklist. Default is FALSE.
#' @param use_curated_blacklist Boolean parameter for using a curated blacklist, default is FALSE.
#' @param verbose For debugging, print out a bunch of possibly useful information.
#' @param invert USE WITH CAUTION! This returns only the variants that would be dropped in the process (opposite of what you want, probably).
#' 
#' @return A MAF format data frame with two new columns indicating the number of occurrences of each variant in the two blacklists.
#' 
#' @import dplyr readr tidyr
#' 
#' @export
#'
#' @examples 
#'
#' #annotate MAF
#' deblacklisted_maf = annotate_ssm_blacklist(grande_maf,
#'                                            seq_type = "genome",
#'                                            genome_build = "hg38")
#'
annotate_ssm_blacklist = function(mutations_df,
                                  seq_type,
                                  tool_name = "slms_3",
                                  tool_version = "1.0",
                                  annotator_name = "vcf2maf",
                                  annotator_version = "1.2",
                                  genome_build = "grch37",
                                  project_base,
                                  blacklist_file_template,
                                  drop_threshold = 4,
                                  return_blacklist = FALSE,
                                  use_curated_blacklist = FALSE,
                                  verbose = FALSE,
                                  invert = FALSE){

  if(missing(seq_type)){
    message("User must specify seq_type of the mutations to select the right blacklist file. More than one seq_type can be specified as a vector if desired.")
    return()
  }

  projection = genome_build

  if(missing(blacklist_file_template)){
    blacklist_template = check_config_value(config::get("resources")$blacklist$template)
  }else{
    blacklist_template = blacklist_file_template
  }

  if(missing(project_base)){
    project_base = check_config_value(config::get("project_base"))
  }
  
  if(!use_curated_blacklist){
    blacklist_files = glue::glue(blacklist_template)
    blacklist_list = list()
    for(b in blacklist_files){
      full_path = paste0(project_base,b)

      #check for missingness
      if(!file.exists(full_path)){
        print(paste("missing: ", full_path))
        message("Cannot find file locally. If working remotely, perhaps you forgot to load your config (see below) or sync your files?")
        message('Sys.setenv(R_CONFIG_ACTIVE = "remote")')
      }

      lifted_blacklist = suppressMessages(readr::read_tsv(full_path, col_names = c("chrpos", "blacklist_count"), col_types = "ci"))

      lifted_blacklist = lifted_blacklist %>%
        separate(chrpos, into = c("Chromosome", "Start_Position"), sep = ":") %>% 
        mutate(Start_Position = as.numeric(Start_Position))
  
      blacklist_list[[b]] = lifted_blacklist
    }

    combined_blacklist = do.call("rbind", blacklist_list)
    
    if(return_blacklist){
      return(combined_blacklist)
    }
  
    #join using chromosome and position
    if(verbose){
      print(head(mutations_df))
      print(head(combined_blacklist))
    }

    if(str_detect(mutations_df$Chromosome, "chr")[1]){
      combined_blacklist = mutate(combined_blacklist, Chromosome = paste0("chr", Chromosome))
    }
    
    mutations_df = left_join(mutations_df, combined_blacklist, by = c("Chromosome", "Start_Position")) %>%
      mutate(blacklist_count = replace_na(blacklist_count, 0))
    dropped = dplyr::filter(mutations_df, blacklist_count > drop_threshold)

    if(verbose){
      if(nrow(dropped) > 0 ){
        ndrop = length(dropped$Tumor_Sample_Barcode)
        message(paste(ndrop, "variants were dropped"))
      }else{
        message("0 variants were dropped")
      }
    }
  
  }else{
    repo_base = check_config_value(config::get("repo_base"))
    full_path = paste0(repo_base, check_config_value(config::get("resources")$curated_blacklist))
    
    additional_blacklist = glue::glue(full_path) %>% 
    readr::read_tsv()

    additional_blacklist = additional_blacklist %>%
      separate(chrpos, into = c("Chromosome", "Start_Position"), sep = ":") %>% 
      mutate(Start_Position = as.numeric(Start_Position))

    mutations_df = left_join(mutations_df, additional_blacklist, by = c("Chromosome", "Start_Position")) %>%
      mutate(blacklist_count = tidyr::replace_na(blacklist_count, 0))

    dropped = dplyr::filter(mutations_df, blacklist_count > drop_threshold)

    if(verbose){
      if(nrow(dropped) > 0 ){
        ndrop = length(dropped$Tumor_Sample_Barcode)
        message(paste(ndrop, "variants were dropped"))
      } else {
        message("0 variants were dropped")
      }
    }
  }
  
  #drop anything that exceeds our threshold but keep NA
  mutations_df = dplyr::filter(mutations_df, is.na(blacklist_count) | blacklist_count < drop_threshold)
  if(invert){
    return(dropped)
  }
  return(mutations_df)
}


#' @title Annotate Recurrent CNVs.
#' 
#' @description Annotates recurrent CNVs from a data frame with CNV data.
#'
#' @details This function takes a data frame with CNVs (`seq_df`) and annotates recurrent CNVs.
#'
#' @param seg_df Data frame of sequences with start and end coordinates.
#' @param seg_file Optional argument to read sequences from file (currently not used in function).
#'
#' @return Nothing.
#'
#' @examples
#' \dontrun{
#' my_segs = get_sample_cn_segments(these_sample_ids = "HTMCP-01-06-00422-01A-01D")
#' annotated = annotate_recurrent_cnv(seg_df = my_segs)
#' }
#' 
annotate_recurrent_cnv = function(seg_df,
                                  seg_file){

  cnv_coord_df = data.frame(chrom = c("18"),
                            start = c(60000000),
                            end = c(61000000),
                            anchor = c("right"),
                            min_cn = c(3),
                            name = "18der_gain")

  cnv_dt = as.data.table(cnv_coord_df)
  seg_dt = as.data.table(seg_df)
  setkey(cnv_dt, chrom, start, end)
  setkey(seg_dt, start, end)
}


#' @title Annotate IGH Breakpoints.
#'
#' @description Add annotation to IGH breakpoints to infer mechanism based on location within IGH.
#'
#' @details Returns a modified bedpe with additional columns.
#' The most useful columns that are added are mechanism (one of CSR, AID, VDJ) and label (NA if unmatched, otherwise one of Emu, Smu, one of the J segments or switch regions).
#'
#' @param annotated_df Previously annotated data frame of SVs.
#' @param genome_build Version of reference build to be used, only grch37 currently accepted.
#'
#' @return A slightly modified bedpe with added columns.
#'
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import dplyr
#' @export
#'
#' @examples
#' manta_sv = get_manta_sv(verbose = FALSE)
#' all_annotated = annotate_sv(sv_data = manta_sv)
#' ig_annotated = annotate_igh_breakpoints(all_annotated)
#'
annotate_igh_breakpoints = function(annotated_df,
                                    genome_build = "grch37"){

  if(genome_build != "grch37"){
    message("Currently, only grch37 is supported")
    return()
  }
  sv_igh = dplyr::filter(annotated_df, partner == "IGH")

  switch_end = 106327893 #hg19 coordinate
  vdj_start = 106330635
  emu_start= 106327853
  emu_end = 106330900
  j_segs = data.frame(start = c(106329408, 106330010, 106330412, 106330789, emu_start, 106323800, 106054077, 106320301, 106304755, 106231780, 106206750, 106172700, 106107740, 106066513, 106089185, 106134348),
                      end = c(106329487, 106330175, 106330579, 106330991, emu_end, 106327600, 106058458, 106322397, 106312100, 106237884, 106209489, 106174976, 106112209, 106068295, 106092889, 106136399),
                      label = c("J6", "J5", "J4", "J3", "Emu", "Smu", "A2", "M", "D", "G3", "G1", "A1", "G2", "E", "G4", "GP"),
                      y = c(-0.5, -0.5, -0.5, -0.5, 0, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5))

  annoying_segments = c("J6", "J3P", "J5", "J4", "J3", "Smu", "M")
  j_segs = mutate(j_segs, start_wide = ifelse(!label %in% annoying_segments, start - 3500, start), end_wide = ifelse(!label %in% annoying_segments, end + 5000, end))
  sv_igh_anno = mutate(sv_igh, mechanism = case_when(chrom2 %in% c("14", "chr14") & start2 < switch_end ~ "CSR",
                                                     chrom2 %in% c("14", "chr14") & start2 < vdj_start ~ "AID",
                                                     chrom2 %in% c("14", "chr14") & start2 >= vdj_start ~ "VDJ",
                                                     chrom1 %in% c("14", "chr14") & start1 < switch_end ~ "CSR",
                                                     chrom1 %in% c("14", "chr14") & start1 < vdj_start ~ "AID",
                                                     chrom1 %in% c("14", "chr14") & start1 >= vdj_start ~ "VDJ",
                                                     TRUE ~ "OTHER"))
  sv_igh_anno = sv_igh_anno %>%
    mutate(IGH_start = ifelse(chrom2 %in% c("14", "chr14"), start2, start1)) %>%
    mutate(IGH_end = IGH_start + 1)

  #assign a higher resolution position by matching to nearest annotation
  sv_igh_anno_dt = as.data.table(sv_igh_anno)

  j_segs_dt = dplyr::filter(j_segs, label != "Emu") %>%
    as.data.table()

  #annotate Emu separately because it overlaps with some annoying J segments
  setkey(sv_igh_anno_dt, IGH_start, IGH_end)
  setkey(j_segs_dt, start_wide, end_wide)
  foverlaps(sv_igh_anno_dt, j_segs_dt, by.x = c("IGH_start", "IGH_end"), by.y = c("start_wide", "end_wide"))
  annotated = foverlaps(sv_igh_anno_dt,j_segs_dt, by.x = c("IGH_start", "IGH_end"), by.y = c("start_wide", "end_wide"))
  annotated[is.na(label) & IGH_start > emu_start & IGH_end < emu_end, "label"] = "Emu"
  return(annotated)
}


#' @title Annotate Driver SSM.
#'
#' @description Retrieve data or use the supplied data frame of genes to get/annotate MAF data indicating which rows are putative driver mutations.
#'
#' @details Provide a maf-like data frame with `maf_df` as the only required parameter. 
#' For information on how to use the additional parameters, refer to the parameter descriptions and function examples. 
#'
#' @param maf_df Data frame of MAF-format mutations.
#' @param lymphoma_type Optional keyword to find genes for annotation (e.g. "BL","DLBCL","FL") is not specifying driver genes.
#' @param driver_genes Optional vector of Hugo_Symbol of genes for coding mutations.
#' @param noncoding_regions Optional named vector of regions to use to further restrict noncoding mutations per gene.
#'
#' @return Incoming MAF, subset to putative driver mutations.
#'
#' @import dplyr
#' @export
#'
#' @examples
#' driver_ssm = annotate_driver_ssm(maf_df = grande_maf,
#'                                  lymphoma_type = "DLBCL",
#'                                  noncoding_regions=c("NFKBIZ"="chr3:101578206-101578365"))
#'
annotate_driver_ssm = function(maf_df,
                               lymphoma_type,
                               driver_genes,
                               noncoding_regions = c("NFKBIZ" = "chr3:101578206-101578365", "HNRNPH1" = "chr5:179,045,946-179,046,683")){

  #get the gene list if missing
  if(missing(driver_genes)){
    driver_genes = GAMBLR.data::lymphoma_genes_lymphoma_genes_v0.0[which(GAMBLR.data::lymphoma_genes_lymphoma_genes_v0.0[[lymphoma_type]] == TRUE),] %>%
      pull(Gene)

    ngene = length(driver_genes)
    message(paste("using", ngene, "genes"))
  }
  if(missing(maf_df)){
    warning("Please specify mutations data frame (MAF-format)")
    return()
  }else{
    #this all assumes the user has left the noncoding variants in their MAF
    kept_ssm = maf_df %>%
      dplyr::filter(Variant_Classification %in% coding_vc) %>%
      dplyr::filter(Hugo_Symbol %in% driver_genes)

    for(gene in names(noncoding_regions)){
      message(paste("adding", unname(noncoding_regions[gene]), "for", gene))
      nc_ssm = maf_df %>%
      dplyr::filter(Variant_Classification == unname(noncoding_regions[gene]))

      if(!is.na(noncoding_regions[gene])){
        #also restrict to coordinates
        chunks = region_to_chunks(noncoding_regions[gene])
        nc_ssm = dplyr::filter(nc_ssm,Start_Position >= chunks$start & Start_Position <= chunks$end)
      }
      kept_ssm = bind_rows(kept_ssm, nc_ssm)
      print(dim(kept_ssm))
    }
  }
  return(kept_ssm)
}


#' @title Annotate SVs.
#'
#' @description Annotate a data frame of SV breakpoints represented in an extended BEDPE format.
#'
#' @details Specify a data frame with SVs (preferably the output from [GAMBLR::get_manta_sv]) to the `sv_df` parameter and get back the same data frame with SV annotations.
#'
#' @param sv_data A data frame of SVs. This should be the output of get_manta_sv. If you aren't using the database backend you can supply your own data frame in the format show below.
#' Most of this data is directly from the bedpe files that are obtained by converting the Manta outputs from VCF.
#' Only the columns for both chromosomes, coordinates and strand plus SOMATIC_SCORE and tumour_sample_id are absolutely required.
#'  CHROM_A  START_A    END_A CHROM_B  START_B    END_B NAME SOMATIC_SCORE STRAND_A STRAND_B TYPE FILTER VAF_tumour VAF_normal DP_tumour DP_normal tumour_sample_id normal_sample_id pair_status
#'   1  1556541  1556547       1  1556664  1556670    .            40        -        -  BND   PASS      0.145          0        55        73  00-14595_tumorA  00-14595_normal     matched
#' @param partner_bed Optional bed-format data frame to use for annotating oncogene partners (e.g. enhancers). required columns are: chrom,start,end,gene
#' @param with_chr_prefix Optionally request that chromosome names are returned with a chr prefix. Default is FALSE.
#' @param collapse_redundant Remove reciprocal events and only return one per event. Default is FALSE.
#' @param return_as Stated format for returned output, default is "bedpe". Other accepted output formats are "bed" and "bedpe_entrez" (to keep entrez_ids for compatibility with portal.R and cBioPortal).
#' @param blacklist A vector of regions to be removed from annotations. Default coordinates are in respect to hg19.
#' @param genome_build Reference genome build parameter, default is grch37.
#'
#' @return A data frame with annotated SVs (gene symbol and entrez ID).
#' 
#' @rawNamespace import(data.table, except = c("last", "first", "between", "transpose"))
#' @import tidyr dplyr stringr
#' @export
#'
#' @examples
#' sv_df = get_manta_sv(verbose = FALSE)
#' annotated_entrez = annotate_sv(sv_data = sv_df,
#'                                with_chr_prefix = FALSE,
#'                                collapse_redundant = FALSE,
#'                                return_as = "bedpe_entrez", 
#'                                genome_build = "grch37")
#'
annotate_sv = function(sv_data,
                       partner_bed,
                       with_chr_prefix = FALSE,
                       collapse_redundant = FALSE,
                       return_as = "bedpe",
                       blacklist = c(60565248, 30303126, 187728894, 101357565, 101359747, 161734970, 69400840, 65217851, 187728889, 187728888,187728892, 187728893,188305164),
                       genome_build = "grch37"){

  bedpe1 = sv_data %>%
    dplyr::select("CHROM_A", "START_A", "END_A", "tumour_sample_id", ends_with("SCORE"), "STRAND_A")

  bedpe2 = sv_data %>%
    dplyr::select("CHROM_B", "START_B", "END_B", "tumour_sample_id", ends_with("SCORE"), "STRAND_B")

  colnames(bedpe1) = c("chrom", "start", "end", "tumour_sample_id", "score", "strand1")
  colnames(bedpe2) = c("chrom", "start", "end", "tumour_sample_id", "score", "strand2")
  suppressWarnings({
    if(any(grepl("chr", bedpe1$chrom))){
      bedpe1 = dplyr::mutate(bedpe1, chrom = str_replace(chrom, "chr", ""))
      bedpe2 = dplyr::mutate(bedpe2, chrom = str_replace(chrom, "chr", ""))
    }
  })
  if(missing(partner_bed)){
    if(genome_build == "hg38"){
      ig_regions = hg38_partners
    }else{
      ig_regions = grch37_partners
    }
  }else{
    ig_regions = partner_bed
    if(!"entrez" %in% colnames(ig_regions)){
      ig_regions$entrez = 0
    }
  }
  if(genome_build == "hg38"){
    oncogene_regions = hg38_oncogene
  }else{
    oncogene_regions = grch37_oncogene
  }
  y = data.table::as.data.table(oncogene_regions)
  data.table::setkey(y, chrom, start, end)

  #use foverlaps to get oncogene SVs
  a = data.table::as.data.table(bedpe1)
  a.onco = data.table::foverlaps(a, y, type = "any", mult = "first") #oncogene-annotated bedpe for the first breakpoints

  b = data.table::as.data.table(bedpe2)
  b.onco = data.table::foverlaps(b, y, type = "any", mult = "first") #oncogene-annotated bedpe for the first breakpoints

  #insist oncogene breakpoints are anchored in an IG or superenhancer region (currently just IG or BCL6)
  #other end of breakpoint
  a.onco.break = a.onco[which(!is.na(a.onco$start)), c("chrom", "i.start", "i.end", "tumour_sample_id", "gene", "entrez", "score", "strand1")]
  b.onco.break = b.onco[which(!is.na(b.onco$start)), c("chrom", "i.start", "i.end", "tumour_sample_id", "gene", "entrez", "score", "strand2")]

  a.partner = b[which(!is.na(a.onco$start)),]
  b.partner = a[which(!is.na(b.onco$start)),]

  y = data.table::as.data.table(ig_regions)
  data.table::setkey(y, chrom, start, end)

  a.ig = data.table::foverlaps(a.partner, y, type = "any", mult = "first")
  b.ig = data.table::foverlaps(b.partner, y, type = "any", mult = "first")

  a.ig = a.ig[,c("chrom", "i.start", "i.end", "strand2", "gene")]
  b.ig = b.ig[,c("chrom", "i.start", "i.end", "strand1", "gene")]

  a.annotated.both = cbind(a.onco.break, a.ig)
  colnames(a.annotated.both) = c("chrom1", "start1", "end1", "tumour_sample_id", "gene", "entrez", "score", "strand1", "chrom2", "start2", "end2", "strand2", "partner")

  b.annotated.both = cbind(b.onco.break, b.ig)
  colnames(b.annotated.both) = c("chrom2", "start2", "end2", "tumour_sample_id", "gene", "entrez", "score", "strand2", "chrom1", "start1", "end1", "strand1", "partner")

  all.annotated = rbind(a.annotated.both, b.annotated.both)

  all.annotated$fusion = dplyr::pull(tidyr::unite(all.annotated, fusion, partner, gene, sep = "-"), fusion)
  all.annotated = dplyr::filter(all.annotated, !fusion %in% c("BCL6-BCL6", "CIITA-CIITA", "FOXP1-FOXP1"))

  all.annotated = dplyr::filter(all.annotated, !start1 %in% blacklist)
  all.annotated = dplyr::filter(all.annotated, !start2 %in% blacklist)

  if(return_as == "bedpe"){
    all.annotated$name = "."
    all.annotated = dplyr::select(all.annotated, chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, tumour_sample_id, gene, partner, fusion)
  }else if(return_as == "bedpe_entrez"){ #necessary for setting up cBioPOrtal instances (setup_study from portal.R)
    all.annotated$name = "."
    all.annotated = dplyr::select(all.annotated, chrom1, start1, end1, chrom2, start2, end2, name, score, strand1, strand2, tumour_sample_id, gene, entrez, partner, fusion)
  }else if(return_as == "bed"){

    #lose the linkage but add a name that somewhat retains it
    if(!grepl("chr", all.annotated$chrom1)){
      all.annotated = all.annotated %>%
        dplyr::mutate(chrom1 = paste0("chr", chrom1))

      all.annotated = all.annotated %>%
        dplyr::mutate(chrom2 = paste0("chr", chrom2))
    }
    bed1 = dplyr::mutate(all.annotated, name = paste(tumour_sample_id, fusion, sep = "_")) %>%
      dplyr::select(chrom1, start1, end1, name, score, strand1)

    bed2 = dplyr::mutate(all.annotated, name = paste(tumour_sample_id, fusion, sep = "_")) %>%
      dplyr::select(chrom2, start2, end2, name, score, strand2)

    colnames(bed1) = c("chrom", "start", "end", "name", "score", "strand")
    colnames(bed2) = c("chrom", "start", "end", "name", "score", "strand")
    return(dplyr::arrange(rbind(bed1, bed2), name))
  }else{
    if(collapse_redundant){
      all.annotated = dplyr::distinct(all.annotated, tumour_sample_id, fusion, .keep_all = TRUE)
    }
  }
  if(with_chr_prefix){

    #add the prefix if necessary
    if(!grepl("chr", all.annotated$chrom1)){
      all.annotated = all.annotated %>%
        dplyr::mutate(chrom1 = paste0("chr", chrom1))

      all.annotated = all.annotated %>%
        dplyr::mutate(chrom2 = paste0("chr", chrom2))
    }
  }
  return(all.annotated)
}


#' @title Annotate mutations target AID motif
#'
#' @description Checks for the presence of mutations at a given AID motif
#'
#' @details In positions that reference allele  "C" has been mutated, two nucleotides before and one nucleotide after it 
#' will be captured. In positions that reference allele  "G" has been mutated, one nucleotide before and two nucleotides after it 
#' will be captured. The other mutations will have "NO" for the sequence. 
#' Then, the forward orientation and reverse compliment of the motif will be generated and their presence in aquired reference
#' allele "C" and reference allele "G" sequence will be checked, respectively. ("NO" for the rows that have different reference
#' alleles) 
#'
#' @param maf MAF data frame (required columns: Reference_Allele, Chromosome, Start_Position, End_Position)
#' @param motif The motif sequence (default is WRCY)
#' @param projection The genome build projection for the variants you are working with (default is grch37)
#' @param fastaPath Can be a path to a FASTA file
#'
#' @return A data frame with two extra columns (seq and motif).
#'
#' @rawNamespace import(IRanges, except = c("merge", "shift", "collapse", "union", "slice", "intersect", "setdiff", "desc", "reduce")) 
#' @rawNamespace import(GenomicRanges, except = c("merge", "shift", "union", "intersect", "setdiff", "reduce"))
#' @import Rsamtools readr dplyr
#' @export
#'
#' @examples
#' annotate_ssm_motif_context (maf = maf,
#'           motif = "WRCY",
#'           projection = "hg38"
#'          )


#This function check that if the motif pattern is in the sequence
annotate_ssm_motif_context <- function(maf,
                     motif = "WRCY",
                     projection = "grch37",
                     fastaPath
                     ){
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
    }
    # It checks for the presence of a local fastaPath
    if (!file.exists(fastaPath)) {
        stop("Failed to find the fasta file")
    }
    # It checks that if the user considers "WRCY" as the motif
    if (!motif == "WRCY"){
        stop("It is only currently available for WRCY. ",
             "To use it for another motif, comment \"if\" code for the motif"
        )
    }
    # Create a reference to an indexed fasta file.
    fasta <- Rsamtools::FaFile(file = fastaPath)
    # This section provides the sequence
    # If the reference allele is C, it will return the one nucleotide after it
    # and 2 nucleotides before it
    sequences <- maf %>%
        dplyr::mutate(
            seq = ifelse(
                Reference_Allele == "C",
                    Rsamtools::getSeq(
                        fasta,
                        GenomicRanges::GRanges(
                            maf$Chromosome,
                            IRanges(
                                start = maf$Start_Position - 2,
                                end = maf$End_Position + 1
                            )
                        )
                    ),
                     # If the reference allele is G, it will return 
                     # 2 nucleotides after it and the one nucleotide before it
                     ifelse(
                         Reference_Allele == "G",
                         Rsamtools::getSeq(
                             fasta,
                             GenomicRanges::GRanges(
                                 maf$Chromosome,
                                 IRanges(
                                     start = maf$Start_Position - 1,
                                     end = maf$End_Position + 2
                                 )
                             )
                          ),
                           # In other cases it returns "NO" for that mutation
                           "NO"
                     )
            )
        )
  
    # This section provides motif and its reverse complement 
    compliment <- c(
        'A'= 'T',
        'T'= 'A',
        'C'= 'G',
        'G'= 'C',
        'U'= 'A',
        'R'= 'Y',
        'Y'= 'R',
        'S'= 'S',
        'W'= 'W',
        'K'= 'M',
        'M'= 'K',
        'B'= 'V',
        'D'= 'H',
        'H'= 'D',
        'V'= 'B',
        'N'= 'N'
    )
    IUPACCodeDict <- c(
        'A'= 'A',  # Adenine
        'C'= 'C',  # Cytosine
        'G'= 'G',  # Guanine
        'T'= 'T',  # Thymine
        'R'= '[AG]',  # A or G
        'Y'= '[CT]',  # C or T
        'S'= '[GC]',  # G or C
        'W'= '[AT]',  # A or T
        'K'= '[GT]',  # G or T
        'M'= '[AC]',  # A or C
        'B'= '[CGT]',  # C or G or T
        'D'= '[AGT]',  # A or G or T
        'H'= '[ACT]',  # A or C or T
        'V'= '[ACG]',  # A or C or G
        'N'= '[ACGT]'  # any base
    )
    word <- motif
    splitWord <- strsplit(word,"")[[1]] # Split the word into its letters
    splitWordLen <- length(splitWord)
    forMotif <- character(splitWordLen) # forMotif, the same length as splitWord
    for (i in seq_along(splitWord)){
        # Convert incomplete nuc specification into their different nucleotides
        if (splitWord[i] %in% names(IUPACCodeDict)){
            forMotif[i] = IUPACCodeDict[splitWord[i]]
        }
    }
    # Collapsing all the letters of forward orientation and make it 
    # into a single string
    strForMotif <- paste(forMotif, collapse = "")
    RevCompMotif <- character(splitWordLen)
    for (i in seq_along(splitWord)){
        if (splitWord[i] %in% names(compliment)){
            # Provide the complement version of the motif
            RevCompMotif[i] = compliment[splitWord[i]]
        }
    }
    # Provide the reverse version of RevCompMotif
    vecRevComp <- rev(RevCompMotif)
    IUPACRevCompMotif <- character(splitWordLen)
    for (i in seq_along(vecRevComp)){
        if (vecRevComp[i] %in% names(IUPACCodeDict)){
            IUPACRevCompMotif[i] = IUPACCodeDict[vecRevComp[i]]
        } 
    }
    # Collapsing all the letters of reverse complement orientation and make it 
    # into a single string
    strRevComp <- paste(IUPACRevCompMotif, collapse = "")
  
    #This section checks for the presence of the motif in the sequence
    # If the reference allele is C, it will check the forward orientation with
    # its sequence
    finder <- sequences %>%
        dplyr::mutate(
            !!motif := ifelse(
              Reference_Allele == "C",
              stringr::str_detect(
                  sequences$seq, strForMotif
              ),
              # If the reference allele is G, it will check the reverse-
              # compliment orientation with its sequence
              ifelse(
                  Reference_Allele == "G",
                  stringr::str_detect(
                      sequences$seq, strRevComp
                  ),
                   # For the rest of mutations it will return "NO"
                   "NO"
              )
            )
        )
    return(finder)
}
