
setwd("~/git/gambl/") # set this path to point to your local clone of the gambl repo
# in your config.yml file set repo_base under the remote config (at the top) to point to this path also
library(GAMBLR)
Sys.setenv(R_CONFIG_ACTIVE= "remote")


config::get("project_base")
#[1] "/Users/rmorin/gambl_results/"


# various tests
# 1. get the metadata for genome and capture samples
all_meta = get_gambl_metadata()
all_meta = mutate(all_meta,pathology=ifelse(pathology=="BLL","BL",pathology))
cap_meta = get_gambl_metadata(seq_type_filter = "capture")

cap_meta = cap_meta %>% dplyr::filter(!sample_id %in% all_meta$sample_id)

#load collated summary results
collated = collate_results(from_cache=TRUE)
collated = collate_results(from_cache=TRUE,seq_type_filter="capture")

#2. get mutations (assuming the files have been synced)

coding_ssm = get_coding_ssm(these_samples_metadata = all_meta,seq_type = "genome")

cap_ssm = get_coding_ssm(these_samples_metadata = cap_meta,seq_type = "capture")

cap_ssm_hg38 = get_coding_ssm(these_samples_metadata = cap_meta,seq_type = "capture",projection = "hg38")


dim(coding_ssm)
dim(cap_ssm)
all_ssm = bind_rows(coding_ssm,cap_ssm)


blacklisted = GAMBLR:::annotate_ssm_blacklist(cap_ssm_hg38,seq_type = "capture",genome_build = "hg38")

#test SV functions
svar_all = get_combined_sv(min_vaf = 0.1,projection = "grch37")
#test annotation
svar_anno = annotate_sv(svar_all)


# test CN functionality

all_seg = get_sample_cn_segments(sample_list = all_meta$sample_id,multiple_samples = T)

# test if we can make a cBioportal instance

cbio_path = paste0(config::get("project_base"),"cbioportal-docker-compose/study/GAMBL_capture_2022/")

samples_included = GAMBLR::setup_study(seq_type = "capture",
                    short_name="GAMBL_capture_2022",
                    human_friendly_name = "GAMBL exomes 2022 edition",
                    project_name="gambl_capture_2022",
                    description = "GAMBL capture data",
                    out_dir = cbio_path)

finalize_study(seq_type="capture",short_name="GAMBL_capture_2022",sample_ids=samples_included,
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
finalize_study(seq_type_filter="genome",short_name="GAMBL_genome_2022",sample_ids=all_samples,
               human_friendly_name = "GAMBL genomes 2022 edition",
               project_name="gambl_genome_2022",
               description = "GAMBL genome data",out_dir = cbio_path,overwrite = TRUE)

session = GAMBLR::get_ssh_session()
test_ssm = get_ssm_by_samples(these_sample_ids = c("14-24534_tumorA","14-24534_tumorB"),
                              ssh_session = session,subset_from_merge = F)


test_ssm1 = get_ssm_by_patients(these_patient_ids =  c("14-24534"),
                              ssh_session = session,subset_from_merge = F)

table(test_ssm$Tumor_Sample_Barcode) #should be the same
table(test_ssm1$Tumor_Sample_Barcode) #should be the same

#copy number functionality

cn_seg = get_sample_cn_segments("14-24534_tumorB",from_flatfile = T)

cn_ssm = assign_cn_to_ssm("14-24534_tumorB",seg_file_source = "battenberg")

pursteenah = estimate_purity(sample_id = "14-24534_tumorB",seg_file_source = "battenberg")
