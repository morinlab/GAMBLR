
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


