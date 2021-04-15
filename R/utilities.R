require(tidyverse)
require(rtracklayer)
require(SRAdb)

#' Title
#'
#' @param sample
#' @param patient
#' @param seq_type
#'
#' @return
#' @export
#'
#' @examples
get_sample_metadata = function(sample,patient,seq_type="tumour"){
  meta = get_gambl_metadata()
  if(seq_type == "tumour"){
    if(!missing(sample)){
      meta %>% filter(sample_id == sample) %>% pull(data_path)
    }
    if(!missing(patient)){
      meta %>% filter(patient_id == patient) %>% pull(data_path)
    }
  }
  return(meta)
}


#' Title
#'
#' @param sample
#' @param patient
#'
#' @return
#' @export
#'
#' @examples
get_bams = function(sample,patient){
  meta = get_gambl_metadata(tissue_status_filter = c("tumour","normal"))
  meta_mrna = get_gambl_metadata(seq_type_filter = "mrna")
  #get all samples for this patient
  if(missing(patient)){
    patient = meta %>% filter(sample_id==sample) %>% pull(patient_id)
  }
  meta_patient = meta %>% filter(patient_id == patient)
  meta_mrna_patient = meta_mrna %>% filter(patient_id == patient)
  build = pull(meta_patient,genome_build) %>% head(1)
  if(build == "hs37d5"){
    igv_build = "hg19"
  }else{
    igv_build = build
  }
  tumour_genome_bams = filter(meta_patient,seq_type == "genome" & tissue_status == "tumour") %>% pull(data_path)
  bam_details = list(igv_build=igv_build, genome_build=build, tumour_bams=tumour_genome_bams)
  normal_genome_bams = filter(meta_patient,seq_type == "genome" & tissue_status == "normal") %>% pull(data_path)
  if(length(normal_genome_bams)){
    bam_details$normal_genome_bams = normal_genome_bams
  }
  if(length(normal_genome_bams)){
    bam_details$normal_genome_bams = normal_genome_bams
  }
  rnaseq_bams = filter(meta_mrna_patient,seq_type == "mrna") %>% pull(data_path)
  if(length(rnaseq_bams)){
    bam_details$rnaseq_bams = rnaseq_bams
  }
  return(bam_details)
}


#' Title
#'
#' @param bams
#' @param genome_build
#' @param region
#' @param padding
#' @param chrom
#' @param start
#' @param end
#' @param sample_id
#' @param out_path
#' @param igv_port
#'
#' @return
#' @export
#'
#' @examples
#' #IMPORTANT: you must be running IGV on the host that is running R and you need to have it listening on a port
#' # The simplest scenario is to run this command on a terminal (if using a Mac), assuming you are using R on gphost10 and you have a ssh config that routes gp10 to that host
#' # ssh -X gp10
#' # then launch IGV (e.e. from a conda installation):
#' # conda activate igv; igv &
#' this_sv = annotated_sv %>% filter(gene=="ETV6")
#' tumour_bam = get_bam_path(sample=this_sv$tumour_sample_id)
#' make_igv_snapshot(chrom=this_sv$chrom2, start=this_sv$start2, end=this_sv$end2, sample_id=this_sv$tumour_sample_id,out_path="~/IGV_snapshots/")
make_igv_snapshot = function(bams,genome_build,region,padding=200,chrom,start,end,sample_id,out_path="/tmp/",igv_port=60506){
  sock= IGVsocket(port = igv_port)
  IGVclear(sock)
  if(missing(region)){
    region=paste0(chrom,":",start-padding,"-",end+padding)
  }
  #region="chr6:90885409-90885410"
  IGVgenome(sock,genome=genome_build)
  IGVgoto(sock,region)
  for(bam_file in bams){
    IGVload(sock,bam_file)
  }
  filename = paste(sample_id,region,"snapshot.png",sep="_")
  IGVsnapshot(sock,fname=filename,dirname=out_path)
  return(paste0(out_path,filename))
}

#' Title
#'
#' @param bedpe_file
#' @param bedpe_df
#' @param target_build
#'
#' @return
#' @export
#'
#' @examples
liftover_bedpe = function(bedpe_file,bedpe_df,target_build="grch37"){
  if(!missing(bedpe_file)){
    original_bedpe = read_tsv(bedpe_file,comment = "##",col_types="cddcddccccccccccccccccc")
  }
  if(!missing(bedpe_df)){
    original_bedpe = bedpe_df
  }
  if(!grepl("chr",original_bedpe$CHROM_A)){
    #add chr prefix
    original_bedpe = original_bedpe %>% mutate(CHROM_A = paste0("chr",CHROM_A)) %>% mutate(CHROM_B = paste0("chr",CHROM_B))
  }
  char_vec = original_bedpe %>% unite(united,sep="\t") %>% pull(united)
  bedpe_obj <- rtracklayer::import(text=char_vec,format="bedpe")
  this_patient = colnames(original_bedpe)[23]
  this_normal = colnames(original_bedpe)[22]
  if(target_build == "grch37" | target_build == "hg19"){
    chain = rtracklayer::import.chain(system.file("extdata","hg38ToHg19.over.chain",package="GAMBLR"))
  }else if(target_build == "grch38" | target_build == "hg38"){
    chain = rtracklayer::import.chain(system.file("extdata","hg19ToHg38.over.chain",package="GAMBLR"))
  }
  colnames(original_bedpe)[1]="CHROM_A"
  original_columns = colnames(original_bedpe)

  first_sv_lifted = rtracklayer::liftOver(bedpe_obj@first,chain)
  second_sv_lifted = rtracklayer::liftOver(bedpe_obj@second,chain)
  no_problem = !((elementNROWS(first_sv_lifted) != 1) | (elementNROWS(second_sv_lifted) != 1))
  first_ok = subset(first_sv_lifted,no_problem)
  second_ok = subset(second_sv_lifted,no_problem)
  first_ok_df = rtracklayer::export(first_ok,format="bed") %>% read_tsv(col_names = c("CHROM_A","START_A","END_A","name_A","score_A","STRAND_A")) %>% select(-score_A) %>% select(-name_A)
  second_ok_df = rtracklayer::export(second_ok,format="bed") %>% read_tsv(col_names = c("CHROM_B","START_B","END_B","name_B","score_B","STRAND_B")) %>% select(-score_B) %>% select(-name_B)
  ok_bedpe = original_bedpe[no_problem,]
  kept_cols = ok_bedpe %>% select(-c("CHROM_A","START_A","END_A","CHROM_B","START_B","END_B","STRAND_A","STRAND_B"))
  fully_lifted = bind_cols(first_ok_df,second_ok_df,kept_cols) %>% select(all_of(original_columns))
  return(fully_lifted)
}

#' Title
#'
#' @param bedpe_paths
#' @param pattern
#' @param out_dir
#'
#' @return
#' @export
#'
#' @examples
read_merge_manta_with_liftover = function(bedpe_paths=c(),pattern="--matched",out_dir){
  to_merge = list()
  for(thispath in bedpe_paths){

    #sample_paths = dir(thispath,pattern=paste0("--",pattern)) #skip unmatched cases for now
    sample_paths = dir(thispath,pattern=pattern) #DEBUGGING
    print(sample_paths)
    #sample_paths = head(sample_paths,15) #for debugging
    for(sp in sample_paths){
      full_path = paste0(thispath,sp,"/somaticSV.bedpe")
      print(paste("working on ", full_path))


      if(grepl("hg38",full_path)){
        print("using liftOver")
        svbed = read_and_liftover_bedpe(full_path) #load and convert to grch37 coordinates
      }else{
        svbed=read_tsv(full_path,comment = "##",col_types="cddcddccccccccccccccccc")
      }

      this_patient = colnames(svbed)[23]
      this_normal = colnames(svbed)[22]
      out_file = paste0(out_dir,"/",this_patient,"--",this_normal,"--hg38Togrch37_sv.tsv")
      print(paste("writing output to",out_file))
      infos = pull(svbed,this_patient)
      infos_n = pull(svbed,this_normal)
      colnames(svbed)[c(1:6)]=c("CHROM_A","START_A","END_A","CHROM_B","START_B","END_B")
      #all_vafs = get.sv.vaf(infos)
      #svbed$VAF = as.numeric(all_vafs)
      svbed$VAF_tumour = sapply(infos,function(x){as.numeric(tail(unlist(strsplit(x,":")),1))})
      svbed$DP_tumour = sapply(infos,function(x){as.numeric(tail(unlist(strsplit(x,":")),2)[1])})
      svbed$VAF_normal = sapply(infos_n,function(x){as.numeric(tail(unlist(strsplit(x,":")),1))})
      svbed$DP_normal = sapply(infos_n,function(x){as.numeric(tail(unlist(strsplit(x,":")),2)[1])})
      svbed$SOMATIC_SCORE = sapply(svbed$INFO_A,function(x){as.numeric(tail(unlist(strsplit(x,"=")),1))})
      #filter on PASS, score, VAF

      #svbed_filt = svbed %>% filter( SCORE > minScore & FILTER == "PASS") %>%
      #  dplyr::select(c(chrom1,start1,end1,chrom2,start2,end2))
      svbed$tumour_sample_id = this_patient
      svbed$normal_sample_id = this_normal
      if(grepl("--unmatched",sp)){
        svbed$pair_status = "unmatched"
      }else{
        svbed$pair_status = "matched"
      }
      print(head(svbed))
      svbed$NAME = "."

      svbed = svbed %>% select(CHROM_A,START_A,END_A,CHROM_B,START_B,END_B,NAME,SOMATIC_SCORE,STRAND_A,STRAND_B,TYPE,FILTER,VAF_tumour,VAF_normal,DP_tumour,DP_normal,tumour_sample_id,normal_sample_id,pair_status)
      #remove chr prefix from both chromosome names
      svbed = svbed %>% mutate(CHROM_A = gsub("chr","",CHROM_A)) %>% mutate(CHROM_B = gsub("chr","",CHROM_B))
      write_tsv(svbed,out_file,col_names=FALSE)
      #to_merge[[this_patient]] = svbed
    }
  }
}

