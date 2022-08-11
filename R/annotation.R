#adding coding_vc to global enviroment
coding_vc = c("Frame_Shift_Del", "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins", "Missense_Mutation", "Nonsense_Mutation", "Nonstop_Mutation", "Splice_Region", "Splice_Site", "Targeted_Region", "Translation_Start_Site")

#' Annotate and auto-drop a MAF data frame with existing blacklists to remove variants that would be dropped during the merge process.
#'
#' @param mutations_df df with mutation data.
#' @param seq_type The seq_type of your mutations if you prefer to apply only the corresponding blacklist (default is to use all available)
#' @param unix_group The unix group (should be the same for all)
#' @param tool_name The tool or pipeline that generated the files (should be the same for all)
#' @param flavour Set to "clustered" if you want to use the blacklist for the new and improved SLMS-3 outputs (otherwise leave empty)
#' @param genome_build The genome build projection for the variants you are working with (default is grch37)
#' @param project_base Optional: A full path to the directory that your blacklist_file_pattern is relative to
#' @param blacklist_file_pattern Optional: A string that contains the relative path to your blacklist file from after the project_base (i.e. results) with any wildcards surrounded with curly braces
#' @param drop_threshold The minimum count from one of the blacklists to drop a variant
#' @param verbose For debugging, print out a bunch of possibly useful information
#' @param invert USE WITH CAUTION! This returns only the variants that would be dropped in the process (opposite of what you want, probably)
#'
#' @return A MAF format data frame with two new columns indicating the number of occurrences of each variant in the two blacklists.
#' @export
#'
#' @examples deblacklisted_maf_df = annotate_ssm_blacklist(original_maf_df)
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
    blacklist_template = config::get("resources")$blacklist$template
  }else{
    blacklist_template = blacklist_file_template
  }
  if(missing(project_base)){
    project_base = config::get("project_base")
  }
  blacklist_files = glue::glue(blacklist_template)
  blacklist_list = list()
  for(b in blacklist_files){
    full_path = paste0(project_base,b)
    lifted_blacklist = read_tsv(full_path, col_names = c("chrpos", "blacklist_count"),col_types="ci")
    lifted_blacklist = lifted_blacklist %>%
      separate(chrpos, into = c("Chromosome", "Start_Position"), sep = ":")

    blacklist_list[[b]] = lifted_blacklist
  }
  combined_blacklist = do.call("rbind", blacklist_list)

  # Collapse variant counts per Start_Position
  combined_blacklist = mutate(combined_blacklist, Start_Position = as.integer(Start_Position)) %>%
    group_by(Start_Position, Chromosome) %>%
    summarize(blacklist_count = sum(blacklist_count)) %>%
    ungroup()

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
  mutations_df = left_join(mutations_df,combined_blacklist,by = c("Chromosome", "Start_Position")) %>%
    mutate(blacklist_count = replace_na(blacklist_count, 0))

  if (use_curated_blacklist){
    repo_base = config::get("repo_base")
    full_path = paste0(repo_base, config::get("results_versioned")$curated_blacklist)
    additional_blacklist = glue(repo_path) %>% read_tsv(c("chrpos", "blacklist_count2"), show_col_types = FALSE)
    additional_blacklist = additional_blacklist %>%
      separate(chrpos, into = c("Chromosome", "Start_Position"), sep = ":")

    additional_blacklist = mutate(additional_blacklist, Start_Position = as.integer(Start_Position)) %>%
      group_by(Start_Position, Chromosome) %>%
      summarize(blacklist_count2 = sum(blacklist_count2)) %>%
      ungroup()

    if(str_detect(mutations_df$Chromosome, "chr")[1]){
      additional_blacklist = mutate(additional_blacklist, Chromosome = paste0("chr", Chromosome))
    }
    mutations_df = left_join(mutations_df,additional_blacklist,by = c("Chromosome", "Start_Position")) %>%
      mutate(blacklist_count2 = replace_na(blacklist_count2, 0))
    dropped = dplyr::filter(mutations_df, blacklist_count2 > drop_threshold)
    mutations_df = dplyr::filter(mutations_df, is.na(blacklist_count2) | blacklist_count2 < drop_threshold)

    if(verbose){
      if(length(dropped) > 0 ){
        ndrop = length(dropped$Tumor_Sample_Barcode)
        message(paste(ndrop, "variants were dropped"))
      } else {
        message("0 variants were dropped")
      }

    }

  }
  dropped = dplyr::filter(mutations_df, blacklist_count > drop_threshold)
  if(verbose){
    if(length(dropped) > 0 ){
      ndrop = length(dropped$Tumor_Sample_Barcode)
      message(paste(ndrop, "variants were dropped"))
    } else {
      message("0 variants were dropped")
    }

  }
  #drop anything that exceeds our threshold but keep NA
  mutations_df = dplyr::filter(mutations_df, is.na(blacklist_count) | blacklist_count < drop_threshold)
  if(invert){
    return(dropped)
  }
  return(mutations_df)
}


#' Annotates recurrent CNVs.
#'
#' @param seq_df Data frame of sequences with start and end coordinates.
#' @param seq_file Optional argument to read sequences from file (currently not used in function).
#'
#' @return Nothing
#' @export
#'
#' @examples
#' annotated = annotate_recurrent_cnv(my_seq_df)
#' annotated = annotate_recurrent_cnv(seq_file = "my_seq_data.txt")
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


#' Add annotation to IGH breakpoints to infer mechanism based on location within IGH.
#'
#' @param annotated_df Previously annotated data frame of SVs.
#' @param genome_build Version of reference build to be used, only grch37 currently accepted.
#'
#' @return A slightly modified bedpe with some added columns. The most useful columns that are added are mechanism (one of CSR, AID, VDJ) and label (NA if unmatched, otherwise one of Emu, Smu, one of the J segments or switch regions)
#' @export
#'
#' @examples
#' all_annotated = get_manta_sv() %>% annotate_sv()
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


#' Retrieve data or use supplied data frame of genes to get/annotate MAF data indicating which rows are putative driver mutations.
#'
#' @param maf_df Data frame of MAF-format mutations (default is to retrieve automatically).
#' @param lymphoma_type Optional keyword to find genes for annotation (e.g. "BL","DLBCL","FL") is not specifying driver genes.
#' @param driver_genes Optional vector of Hugo_Symbol of genes for coding mutations.
#' @param noncoding_regions Optional named vector of regions to use to further restrict noncoding mutations per gene.
#'
#' @return Driver ssm kepts from maf.
#' @export
#'
#' @examples
#' driver_ssm = annotate_driver_ssm(maf_df = maf, lymphoma_type = "DLBCL", noncoding_regions=c("NFKBIZ"="chr3:101578206-101578365"))
#'
annotate_driver_ssm = function(maf_df,
                               lymphoma_type,
                               driver_genes,
                               noncoding_regions = c("NFKBIZ" = "chr3:101578206-101578365", "HNRNPH1" = "chr5:179,045,946-179,046,683")){

  #get the gene list if missing
  if(missing(driver_genes)){
    driver_genes = lymphoma_genes[which(lymphoma_genes[[lymphoma_type]] == TRUE),] %>%
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

    for(gene in names(include_noncoding)){
      message(paste("adding", unname(include_noncoding[gene]), "for", gene))
      nc_ssm = maf_df %>%
      dplyr::filter(Variant_Classification == unname(include_noncoding[gene]))

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


#' Annotate a data frame of SV breakpoints represented in an extended BEDPE format.
#'
#' @param sv_data A data frame of SVs. This should be the output of get_manta_sv. If you aren't using the database backend you can supply your own data frame in the format show below.
#' Most of this data is directly from the bedpe files that are obtained by converting the Manta outputs from VCF.
#' Only the columns for both chromosomes, coordinates and strand plus SOMATIC_SCORE and tumour_sample_id are absolutely required.
#'  CHROM_A  START_A    END_A CHROM_B  START_B    END_B NAME SOMATIC_SCORE STRAND_A STRAND_B TYPE FILTER VAF_tumour VAF_normal DP_tumour DP_normal tumour_sample_id normal_sample_id pair_status
#'   1  1556541  1556547       1  1556664  1556670    .            40        -        -  BND   PASS      0.145          0        55        73  00-14595_tumorA  00-14595_normal     matched
#' @param partner_bed Optional bed-format data frame to use for annotating oncogene partners (e.g. enhancers). required columns are: chrom,start,end,gene
#' @param with_chr_prefix Optionally request that chromosome names are returned with a chr prefix.
#' @param collapse_redundant Remove reciprocal events and only return one per event.
#' @param return_as Stated format for returned output, default is bedpe. Other accepted output format is bed and bedpe_entrez (to keep entrez_ids for compatabillity with portal.R and cBioPortal).
#' @param blacklist A list of regions to be removed from annotations. Default coordinates are in respect to hg19.
#' @param genome_build Reference genome build parameter, default is grch37 (hg38 is also accepted).
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
#' annotated_entrez = annotate_sv(sv_data = sv_df, with_chr_prefix = FALSE, collapse_redundant = FALSE, return_as = "bedpe_entrez", genome_build = "grch37")
#'
annotate_sv = function(sv_data,
                       partner_bed,
                       with_chr_prefix = FALSE,
                       collapse_redundant = FALSE,
                       return_as = "bedpe_entrez",
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
