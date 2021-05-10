#' Generate a colourful multi-panel overview of hypermutation in regions of interest across many samples
#'
#' @param regions_bed Bed file with chromosome coordinates, should contain columns chr, start, end, name (with these exact names)
#' @param regions_to_display Optional vector of names from default regions_bed to use
#' @param metadata A metadata file already subsetted and arranged on the order you want the samples vertically displayed
#' @param classification_column optional. Override default column for assigning the labels used for colouring in the figure.
#'
#' @return nothing
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
ashm_multi_rainbow_plot = function(regions_bed,regions_to_display,exclude_classifications,metadata,custom_colours,classification_column="lymphgen",maf_data){
  table_name = config::get("results_tables")$ssm
  db=config::get("database_name")
  #get the mutations for each region and combine
  #regions_bed should contain chr, start, end, name (with these exact names)
  if(missing(metadata)){
    metadata = get_gambl_metadata()
    meta_arranged = arrange(metadata,pathology_rank,lymphgen)
  }else{
    meta_arranged = metadata #assume the user already arranged it the way they wanted
  }
  if(!missing(exclude_classifications)){
    meta_arranged = filter(meta_arranged,!get(classification_column) %in% exclude_classifications)
  }
  if(missing(regions_bed)){
    regions_bed= grch37_ashm_regions
    regions_bed = mutate(regions_bed,regions=paste0(chr_name,":",hg19_start,"-",hg19_end))
    regions_bed = mutate(regions_bed,name=paste0(gene,"-",region))
  }else{
    regions_bed = mutate(regions_bed,regions=paste0(chr,":",start,"-",end))
  }

  names=pull(regions_bed,name)
  names = c(names,"NFKBIZ-UTR","MAF","PAX5","WHSC1","CCND1",
            "FOXP1-TSS1","FOXP1-TSS2","FOXP1-TSS3","FOXP1-TSS4","FOXP1-TSS5","BCL6","IGH","IGL","IGK","PVT1","BCL2") #add some additional regions of interest
  regions = pull(regions_bed,regions)
  regions = c(regions,"chr3:101578214-101578365","chr16:79627745-79634622","chr9:36898851-37448583","chr4:1867076-1977887","chr11:69451233-69460334","chr3:71623481-71641671","chr3:71532613-71559445","chr3:71343345-71363145","chr3:71167050-71193679","chr3:71105715-71118362",
              "chr3:187406804-188522799","chr14:106144562-106344765","chr22:23217074-23250428","chr2:89073691-89320640",
              "chr8:128774985-128876311","chr18:60982124-60990180")
  regions_bed = data.frame(regions=regions,names=names)
  regions_bed = filter(regions_bed,names %in% regions_to_display)
  regions = pull(regions_bed,regions)
  names=pull(regions_bed,names)
  if(missing(maf_data)){
    region_mafs = lapply(regions,function(x){get_ssm_by_region(region=x,streamlined = TRUE)})
  }else{
    region_mafs = lapply(regions,function(x){get_ssm_by_region(region=x,streamlined=TRUE,maf_data=maf_data)})
  }
  tibbled_data = tibble(region_mafs, region_name = names)
  unnested_df = tibbled_data %>% unnest_longer(region_mafs)
  unlisted_df = mutate(unnested_df,start=region_mafs$Start_Position,sample_id=region_mafs$Tumor_Sample_Barcode) %>%
    select(start,sample_id,region_name)



  meta_arranged$classification = factor(meta_arranged[,classification_column],levels=unique(meta_arranged[,classification_column]))
  muts_anno = left_join(unlisted_df,meta_arranged)
  muts_first =  select(muts_anno,start,region_name) %>% group_by(region_name) %>% arrange(start) %>% filter(row_number()==1)
  eg = expand_grid(start=pull(muts_first,start),sample_id=pull(meta_arranged,sample_id))
  eg = left_join(eg,muts_first)

  #concatenate expanded frame of points with original mutation data
  real_and_fake = bind_rows(unlisted_df,eg)
  muts_anno = left_join(real_and_fake,meta_arranged) #%>% filter(!is.na(get(classification_column)))

  muts_anno$sample_id= factor(muts_anno$sample_id,levels=meta_arranged$sample_id)

  if(!missing(regions_to_display)){
    muts_anno = filter(muts_anno,region_name %in% regions_to_display)
  }
  #make the plot
  if(missing(custom_colours)){
    p = muts_anno %>%
      ggplot() + geom_point(aes(x=start,y=sample_id,colour=classification),alpha=0.4,size=0.6) +
      theme(axis.text.y=element_blank()) +
      facet_wrap(~region_name,scales="free_x") + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))

    print(p)
  }else{
    #testing manual colouring
    p = muts_anno %>%
      ggplot() + geom_point(aes(x=start,y=sample_id,colour=classification),alpha=0.4,size=0.6) +
      theme(axis.text.y=element_blank()) + scale_colour_manual(values=custom_colours) +
      facet_wrap(~region_name,scales="free_x") + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
    print(p)
  }
}

#' Create a genome-wide copy number plot for one sample and (optionally) display mutation VAF
#'
#' @param this_sample
#' @param just_segments
#' @param genes_to_label optional. Provide a list of genes to label (if mutated). Can only be used with coding_only (see below)
#' @param coding_only optional. Set to TRUE to restrict to plotting only coding mutations
#' @param from_flatfile
#'
#' @return nothing
#' @export
#' @import tidyverse DBI RMariaDB
#'
#' @examples
copy_number_vaf_plot = function(this_sample,just_segments=FALSE,coding_only=FALSE,one_chrom,
                                genes_to_label,from_flatfile=FALSE,use_augmented_maf=FALSE){
  chrom_order=factor(c(1:22,"X"))
  cn_colours = get_gambl_colours(classification = "copy_number")
  maf_and_seg = assign_cn_to_ssm(this_sample=this_sample,coding_only=coding_only,from_flatfile=from_flatfile,use_augmented_maf=use_augmented_maf)
  vaf_cn_maf = maf_and_seg[["maf"]]
  vaf_cn_maf = mutate(vaf_cn_maf,CN=case_when(LOH == "1" & CN == 2 ~ "nLOH",
                                              TRUE ~ as.character(CN)))
  if(!missing(one_chrom)){
    vaf_cn_maf = filter(vaf_cn_maf,Chromosome == one_chrom)
  }
  #vaf_cn_maf = mutate(vaf_cn_maf,CN=as.character(CN))
  #use_neut_loh=TRUE

  if(just_segments){
    #I realized this is ugly
    cn_seg = maf_and_seg[["seg"]]
    cn_seg = mutate(cn_seg,CN_segment = as.numeric(CN),CN = as.character(CN))
    ggplot(cn_seg) +
      geom_segment(data=cn_seg,aes(x=Start_Position,xend=End_Position,y=CN_segment,yend=CN_segment)) +
      facet_wrap(~factor(Chromosome,levels=chrom_order),scales="free_x") +
      scale_colour_manual(values = cn_colours) +
      theme_minimal() + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
  }else{
    if(coding_only){
      if(missing(genes_to_label)){
        p = mutate(vaf_cn_maf,vaf=t_alt_count/(t_ref_count+t_alt_count)) %>% ggplot() +
          geom_point(aes(x=Start_Position,y=vaf,colour=CN),alpha=0.6,size=2) +
          scale_colour_manual(values = cn_colours) +
          facet_wrap(~factor(Chromosome,levels=chrom_order),scales="free_x") +
          theme_minimal() + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
        p + ggtitle(this_sample)
      }else{
        #label any mutations that intersect with our gene list
        plot_genes = vaf_cn_maf %>% filter(Hugo_Symbol %in% my_genes)

        p = mutate(vaf_cn_maf,vaf=t_alt_count/(t_ref_count+t_alt_count)) %>% ggplot() +
          geom_point(aes(x=Start_Position,y=vaf,colour=CN),size=2) +
          geom_text(data=plot_genes,aes(x=Start_Position,y=0.8,label=Hugo_Symbol),size=3,angle=90) +
          scale_colour_manual(values = cn_colours) +
          facet_wrap(~factor(Chromosome,levels=chrom_order),scales="free_x") + ylim(c(0,1)) +
          theme_minimal() + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
        p + ggtitle(this_sample)
      }
    }else{

      p = mutate(vaf_cn_maf,vaf=t_alt_count/(t_ref_count+t_alt_count)) %>% ggplot() +
        geom_point(aes(x=Start_Position,y=vaf,colour=CN),alpha=0.6,size=0.2) +
        scale_colour_manual(values = cn_colours) +
        facet_wrap(~factor(Chromosome,levels=chrom_order),scales="free_x") +
        theme_minimal() + guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
      p + ggtitle(this_sample)
    }
  }
}


#TODO migrate viz/plotting functions that don't directly rely on the database to a separate file
#' Make a rainbow plot of all mutations in a region, ordered and coloured by metadata
#'
#' @param mutations_maf A data frame containing mutations (MAF format) within a region of interest (i.e. use the get_ssm_by_region)
#' @param metadata should be a data frame with sample_id as a column that should match Tumor_Sample_Barcode in the database
#' @param classification_column The name of the metadata column to use for ordering and colouring samples
#' @param bed Optional data frame specifying the regions to annotate (required columns: start, end, name)
#'
#' @return ggplot2 object
#' @export
#' @import tidyverse DBI RMariaDB dbplyr
#'
#' @examples
#' # basic usage
#' ashm_rainbow_plot(mutations_maf=get_ssm_by_region(region="chr8:128806578-128806992"),metadata=get_gambl_metadata())
#' # advanced usages
#' mybed = data.frame(start=c(128806578,128805652,128748315), end=c(128806992,128809822,128748880), name=c("TSS","enhancer","MYC-e1"))
#' ashm_rainbow_plot(mutations_maf=my_mutations,metadata=my_metadata,bed=mybed)
ashm_rainbow_plot = function(mutations_maf,
                             metadata,exclude_classifications,drop_unmutated=FALSE,
                             classification_column,
                             bed,region,custom_colours,hide_ids=TRUE){
  table_name = config::get("results_tables")$ssm
  db=config::get("database_name")
  if(!missing(region)){

    region = gsub(",","",region)
    split_chunks = unlist(strsplit(region,":"))
    chromosome = split_chunks[1]

    startend = unlist(strsplit(split_chunks[2],"-"))
    qstart=as.numeric(startend[1])
    qend=as.numeric(startend[2])

    if(missing(mutations_maf)){
      mutations_maf = get_ssm_by_region(region=region,streamlined = TRUE)
    }else{
      #ensure it only contains mutations in the region specified
      mutations_maf = get_ssm_by_region(region=region,streamlined = TRUE,maf_data = mutations_maf)
    }

  }
  if(!missing(classification_column)){
    meta_arranged = arrange(metadata,pathology_rank,lymphgen) #fix this to use the actual column
    if(!missing(exclude_classifications)){
      meta_arranged = filter(meta_arranged,!get(classification_column) %in% exclude_classifications)
    }
  }else{
    classification_column = "lymphgen"
    meta_arranged = metadata
  }


  mutation_positions = mutations_maf %>%
    dplyr::select(Tumor_Sample_Barcode,Start_Position) %>% as.data.frame()
  mutated_cases = pull(mutation_positions,Tumor_Sample_Barcode) %>% unique()
  if(drop_unmutated){
    meta_arranged = meta_arranged %>% filter(sample_id %in% mutated_cases)
  }
  #add a fake mutation at the start position for each sample to ensure every sample shows up
  fake_mutations = data.frame(Tumor_Sample_Barcode=pull(metadata,sample_id),Start_Position=qstart-1000)
  mutation_positions = rbind(mutation_positions,fake_mutations)

  meta_arranged$classification = factor(meta_arranged[,classification_column],levels=unique(meta_arranged[,classification_column]))

  muts_anno = dplyr::left_join(mutation_positions,meta_arranged,by=c("Tumor_Sample_Barcode" = "sample_id"))


  muts_anno$sample_id = factor(muts_anno$Tumor_Sample_Barcode,levels=unique(meta_arranged$sample_id))
  #

  if(missing(custom_colours)){
    p = ggplot(muts_anno) +
      geom_point(aes(x=Start_Position,y=sample_id,colour=classification),alpha=0.4)
  }else{
    p = ggplot(muts_anno) + geom_point(aes(x=Start_Position,y=sample_id,colour=classification),alpha=0.4) +
      scale_colour_manual(values=custom_colours)
  }
  if(missing(bed)){
    p +
      guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))
  }else{
    bed = bed %>% mutate(size = end - start) %>% mutate(midpoint = start + size/2)
    height = length(unique(meta_arranged$sample_id)) + 8
    p = p + geom_rect(data=bed, aes(xmin = start, xmax = end, ymin = 0, ymax = height+5),alpha=0.1) +
      geom_text(data=bed,aes(x = midpoint , y= height ,label=name),size = 2.5,angle=90) +
      guides(color = guide_legend(reverse = TRUE,override.aes = list(size = 3)))

  }
  if(hide_ids){
    p= p + theme(axis.text.y=element_blank())
  }else{
    p = p + theme(axis.text.y = element_text(size = 5))
  }
  return(p)
}

#' Title
#'
#' @param mafs
#' @param sample_id
#' @param genes
#' @param show_noncoding
#' @param detail
#'
#' @return
#' @export
#' @import tidyverse
#'
#' @examples
plot_multi_timepoint = function(mafs,sample_id,genes,show_noncoding=FALSE,detail){
  tp = c("A","B","C")
  title = paste(sample_id,detail,sep="\n")
  i = 1
  for (i in c(1:length(mafs))){
    maf_file = mafs[i]
    time_point=tp[i]
    print(paste(maf_file,time_point))
  }
  if(length(mafs)==2){
    A.maf = fread_maf(mafs[1])
    B.maf = fread_maf(mafs[2])

    A.maf = A.maf %>% dplyr::select(c(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Tumor_Sample_Barcode,HGVSp_Short,t_ref_count,t_alt_count))  %>%
      mutate(VAF=t_alt_count/(t_ref_count+t_alt_count)) %>% mutate(time_point=1) %>% mutate(coord=paste(Chromosome,Start_Position,sep=":"))
    B.maf = B.maf %>% dplyr::select(c(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Tumor_Sample_Barcode,HGVSp_Short,t_ref_count,t_alt_count))  %>%
      mutate(VAF=t_alt_count/(t_ref_count+t_alt_count)) %>% mutate(time_point=2) %>% mutate(coord=paste(Chromosome,Start_Position,sep=":"))
    all.maf=rbind(A.maf,B.maf)
    if(show_noncoding){
      coding.maf = subset_regions(all.maf,shm_regions)
    }
    else{
      coding.maf = filter(all.maf,!Variant_Classification %in% c("Silent","RNA","IGR","Intron","5'Flank","3'Flank","5'UTR"))
      #keep certain 3' UTR mutations, toss the rest
      coding.maf = filter(coding.maf,(Variant_Classification == "3'UTR" & Hugo_Symbol == "NFKBIZ") | (Variant_Classification != "3'UTR"))
    }
    A.rows = which(coding.maf$time_point==1)
    A.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==1 & VAF == 0), "coord"]))
    A.zero = filter(coding.maf,coord %in% A.zero.coords)

    B.rows = which(coding.maf$time_point==2)
    B.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==2 & VAF == 0), "coord"]))
    B.zero = filter(coding.maf,coord %in% B.zero.coords)
    coding.maf$category = "trunk"
    coding.maf[which(coord %in% A.zero.coords),"category"]="not-A"
    coding.maf[which(coord %in% B.zero.coords),"category"]="not-B"



    #actually this is changed in eitehr direction, not just gained
    just_gained_lg_all = filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" )
    #just_gained_lg = filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" & time_point !=2) %>%
    #  mutate(time_point = time_point +0.4)

    just_gained_lg = filter(coding.maf, Hugo_Symbol %in% genes & category != "trunk" & time_point == 2 & VAF > 0) %>%
      mutate(time_point = time_point +0.4)
    print(just_gained_lg)
    just_trunk = filter(coding.maf, Hugo_Symbol %in% genes & category == "trunk" & time_point ==1) %>%
      mutate(time_point = time_point-0.4)
    ggplot(coding.maf,aes(x=time_point,y=VAF,group=coord,colour=category)) +
      geom_point() + geom_line(alpha=0.5) +
      geom_text_repel(data=just_gained_lg,aes(label=Hugo_Symbol),size=4,segment.linetype=0) +
      geom_text_repel(data=just_trunk,aes(label=Hugo_Symbol),size=4,segment.linetype=0) +
      ggtitle(title) +
      theme_minimal()
    return(TRUE)
  }
  if(length(mafs)==3){
    A.maf = fread_maf(mafs[1])
    B.maf = fread_maf(mafs[2])
    C.maf = fread_maf(mafs[3])
    A.maf = A.maf %>% dplyr::select(c(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Tumor_Sample_Barcode,HGVSp_Short,t_ref_count,t_alt_count))  %>%
      mutate(VAF=t_alt_count/(t_ref_count+t_alt_count)) %>% mutate(time_point=1) %>% mutate(coord=paste(Chromosome,Start_Position,sep=":"))
    B.maf = B.maf %>% dplyr::select(c(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Tumor_Sample_Barcode,HGVSp_Short,t_ref_count,t_alt_count))  %>%
      mutate(VAF=t_alt_count/(t_ref_count+t_alt_count)) %>% mutate(time_point=2) %>% mutate(coord=paste(Chromosome,Start_Position,sep=":"))
    C.maf = C.maf %>% dplyr::select(c(Hugo_Symbol,Chromosome,Start_Position,End_Position,Variant_Classification,Tumor_Sample_Barcode,HGVSp_Short,t_ref_count,t_alt_count))  %>%
      mutate(VAF=t_alt_count/(t_ref_count+t_alt_count)) %>% mutate(time_point=3) %>% mutate(coord=paste(Chromosome,Start_Position,sep=":"))
    all.maf=rbind(A.maf,B.maf,C.maf)
    if(show_noncoding){
      coding.maf = subset_regions(all.maf,shm_regions)
    }
    else{
      coding.maf = filter(all.maf,!Variant_Classification %in% c("Silent","RNA","IGR","Intron","5'Flank","3'Flank","5'UTR"))
      #keep certain 3' UTR mutations, toss the rest
      coding.maf = filter(coding.maf,(Variant_Classification == "3'UTR" & Hugo_Symbol == "NFKBIZ") | (Variant_Classification != "3'UTR"))
    }

    A.rows = which(coding.maf$time_point==1)
    A.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==1 & VAF == 0), "coord"]))
    A.zero = filter(coding.maf,coord %in% A.zero.coords)

    B.rows = which(coding.maf$time_point==2)
    B.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==2 & VAF == 0), "coord"]))
    B.zero = filter(coding.maf,coord %in% B.zero.coords)

    C.rows = which(coding.maf$time_point==3)
    C.zero.coords = pull(unique(coding.maf[which(coding.maf$time_point==3 & VAF == 0), "coord"]))
    C.zero = filter(coding.maf,coord %in% C.zero.coords)

    coding.maf$category = "trunk"
    coding.maf[which(coord %in% A.zero.coords),"category"]="not-A"
    coding.maf[which(coord %in% B.zero.coords),"category"]="not-B"
    coding.maf[which(coord %in% C.zero.coords),"category"]="not-C"

    just_gained_lg_all = filter(coding.maf, Hugo_Symbol %in% lg & category != "trunk" )
    just_gained_lg = filter(coding.maf, Hugo_Symbol %in% lg & category != "trunk" & time_point ==3) %>%
      mutate(time_point = time_point +0.4)

    just_trunk = filter(coding.maf, Hugo_Symbol %in% lg & category == "trunk" & time_point ==1) %>%
      mutate(time_point = time_point-0.4)
    ggplot(coding.maf,aes(x=time_point,y=VAF,group=coord,colour=category)) +
      geom_point() + geom_line(alpha=0.3) +
      geom_text_repel(data=just_gained_lg,aes(label=Hugo_Symbol),size=4,segment.linetype=0) +
      geom_text_repel(data=just_trunk,aes(label=Hugo_Symbol),size=4,segment.linetype=0) +
      ggtitle(title) +
      theme_minimal()
    return(TRUE)
  }
  return(FALSE)
}
