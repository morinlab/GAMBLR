
# in your config.yml file set repo_base under the remote config (at the top) to point to this path also
library(GAMBLR)
library(tidyverse)
#library(ssh)
#skip this next line if you are running on the GSC network
Sys.setenv(R_CONFIG_ACTIVE= "remote") 
#Even when working remotely, this step should be handled by the bundled Rprofile but users are responsible for copying that to the right location

check_gamblr_config() #check for missing local files. If a file is missing you need to figure out why. Most common explanations:
# 1. you haven't yet run get_gambl_results.smk and need to do that
# 2. someone added a new file to GAMBLR but hasn't updated the snakefile to support the syncing of that file locally (or telling this function it's not expected). 
# To avoid future confusion this should be addressed by updating either this function, the function it calls or the snakefile. 
# 3. 

session = get_ssh_session() #only run this if you have an active VPN connection
check_gamblr_config(compare_timestamps = T,ssh_session=session)


config::get("project_base")
# should point to an existing path on your local computer where you want the data to be stored
#[1] "/Users/rmorin/gambl_results/"
config::get("repo_base")
# must refer to the location of the local clone of the gambl repo (morinlab/gambl)
#[1] "/Users/rmorin/git/gambl/"

# various tests
# 1. get the metadata for genome and capture samples
all_meta = get_gambl_metadata()
all_meta = mutate(all_meta,pathology=ifelse(pathology=="BLL","BL",pathology))
genome_meta = get_gambl_metadata(seq_type_filter = c("genome","capture"),seq_type_priority = "genome") %>% dplyr::filter(seq_type=="genome")
cap_meta = get_gambl_metadata(seq_type_filter = c("genome","capture"),seq_type_priority = "genome") %>% dplyr::filter(seq_type=="capture")

dplyr::filter(genome_meta,Tumor_Sample_Barcode %in% cap_meta$Tumor_Sample_Barcode) %>% nrow()



#load collated summary results
collated = collate_results(from_cache=TRUE)
collated = collate_results(from_cache=TRUE,seq_type_filter="capture")
dim(collated)
# this currently  prints out:
# 2058 27

#2. get mutations (assuming the files have been synced)
# all these should proceed without error
coding_ssm = get_coding_ssm(these_samples_metadata = genome_meta,seq_type = "genome")

cap_ssm = get_coding_ssm(these_samples_metadata = cap_meta,seq_type = "capture")

cap_ssm_hg38 = get_coding_ssm(these_samples_metadata = cap_meta,seq_type = "capture",projection = "hg38")


dim(coding_ssm)
# [1] 179807     45
dim(cap_ssm)
# [1] 626822     45

coding_ssm$Tumor_Sample_Barcode[which(coding_ssm$Tumor_Sample_Barcode %in% cap_ssm$Tumor_Sample_Barcode)]
# no overlap of IDs between the two because we used seq_type_priority

all_ssm = bind_rows(coding_ssm,cap_ssm)


blacklisted = GAMBLR:::annotate_ssm_blacklist(cap_ssm_hg38,seq_type = "capture",genome_build = "hg38")

grande_gambl_ssm = dplyr::filter(coding_ssm,Tumor_Sample_Barcode %in% grande_maf$Tumor_Sample_Barcode)

dim(grande_gambl_ssm)
#[1] 5865   45
dim(grande_maf)
#[1] 12251   125

#test SV functions
svar_all = get_combined_sv(min_vaf = 0.1,projection = "grch37")

dim(svar_all)
# [1] 450641     19

#test annotation
svar_anno = annotate_sv(svar_all)


# test CN functionality

all_seg = get_sample_cn_segments(sample_list = all_meta$sample_id,multiple_samples = T)
#dim(all_seg)
#[1] 376588      7

# Test remote functionality over ssh

session = GAMBLR::get_ssh_session() # If your local machine username doesn't
# match your GSC username, you need to add "<your_username>@gphost01.bcgsc.ca"
# as the only argument to this function.
test_ssm = get_ssm_by_samples(these_sample_ids = c("14-24534_tumorA","14-24534_tumorB"),
                              subset_from_merge = T)

table(test_ssm$Tumor_Sample_Barcode)

#14-24534_tumorA 14-24534_tumorB 
#11097           12904 

test_ssm1 = get_ssm_by_patients(these_patient_ids =  c("14-24534"),
                                subset_from_merge = F)

table(test_ssm1$Tumor_Sample_Barcode) #should be the same

#14-24534_tumorA 14-24534_tumorB 
#11097           12904 


#copy number functionality

cn_seg = get_sample_cn_segments("14-24534_tumorB",from_flatfile = T)

#should automagically get the single seg file and MAF for the sample
cn_ssm = assign_cn_to_ssm("14-24534_tumorB", seg_file_source = "battenberg")


pursteenah = estimate_purity(this_sample_id = "14-24534_tumorB", seg_file_source = "battenberg")

pursteenah$sample_purity_estimation


#pure=estimate_purity(this_sample_id="12-17272_tumorB",seg_file_source = "battenberg")
#fetching: slms-3                                                                                                                                   
#using flatfile: /projects/nhl_meta_analysis_scratch/gambl/results_local/gambl/battenberg_current/99-outputs/seg/genome--grch37/12-17272_tumorB--12-17272_normal_subclones.igv.seg
#Rows: 45 Columns: 6                                                                                                                                
#── Column specification ────────────────────────────────────────────────────────────────────────────────────────
#Delimiter: "\t"
#chr (2): ID, chrom
#dbl (4): start, end, LOH_flag, log.ratio
#ℹ Use `spec()` to retrieve the full column specification for this data.
#ℹ Specify the column types or set `show_col_types = FALSE` to quiet this message.
#> pure$sample_purity_estimation
#[1] 0.3882176

# advanced functions


coding_ssm_status = get_coding_ssm_status(gene_symbols = c("EZH2","CREBBP","KMT2D"),these_samples_metadata = cap_meta,maf_data = cap_ssm)



# test if we can make a cBioportal instance

cbio_path = paste0(config::get("project_base"),"cbioportal-docker-compose/study/GAMBL_capture_2022/")

samples_included = GAMBLR::setup_study(seq_type_filter = "capture",
                    short_name="GAMBL_capture_2022",
                    human_friendly_name = "GAMBL exomes 2022 edition",
                    project_name="gambl_capture_2022",
                    description = "GAMBL capture data",
                    out_dir = cbio_path)

finalize_study(seq_type_filter="capture",short_name="GAMBL_capture_2022",sample_ids=samples_included,
               human_friendly_name = "GAMBL exomes 2022 edition",
               project_name="gambl_capture_2022",
               description = "GAMBL capture data",out_dir = cbio_path,overwrite = TRUE)



cbio_path = paste0(config::get("project_base"),"cbioportal-docker-compose/study/GAMBL_genome_2022/")

samples_included = GAMBLR::setup_study(seq_type_filter = "genome",
                                       short_name="GAMBL_genome_2022",
                                       human_friendly_name = "GAMBL genomes 2022 edition",
                                       project_name="gambl_genome_2022",
                                       description = "GAMBL genome data",
                                       out_dir = cbio_path)


fusion_samples = setup_fusions(human_friendly_name = "GAMBL genomes 2022 edition",
              project_name = "gambl_genome_2022",
              description = "GAMBL genome data",
              out_dir=cbio_path)
  
all_samples = unique(c(samples_included,fusion_samples))
finalize_study(seq_type_filter="genome", short_name="GAMBL_genome_2022", these_sample_ids=all_samples,
               human_friendly_name = "GAMBL genomes 2022 edition",
               project_name="gambl_genome_2022",
               description = "GAMBL genome data", out_dir = cbio_path, overwrite = TRUE)

