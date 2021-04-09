## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup,warning=FALSE,message=FALSE----------------------------------------
library(GAMBLR)
require(tidyverse)
require(maftools)

## ----metadata-----------------------------------------------------------------

my_meta = get_gambl_metadata()

#reduce to some of the more useful metadata fields

my_meta = my_meta %>% select(sample_id,biopsy_id,myc_ba,cohort,pathology)

print(head(my_meta))


## ----manta--------------------------------------------------------------------

myc_locus_sv = get_manta_sv(region="8:128723128-128774067",pass=FALSE) 

#we can override default that requires SV to meet the Manta "Pass" filtering criterion

print(head(myc_locus_sv))


## ----get_seg_data-------------------------------------------------------------

my_segments = get_cn_segments(chromosome="4",qstart=83274467,qend=83295149)
print(head(my_segments))

deleted_segments = my_segments %>% filter(log.ratio<0) #use some lower value if you want to be more stringent

annotated_segments = left_join(deleted_segments,my_meta,by=c("ID" = "sample_id"))

annotated_segments %>% pull(cohort) %>% table()



## ----ssms_and_maftools--------------------------------------------------------

all_ssms = get_ssm_by_gene(gene_symbol=c("CCND3"),coding_only = TRUE)

all_ssms = all_ssms %>% as.data.frame()

#make a MAFtools object and plot a lollipop plot

maf_obj = read.maf(all_ssms)
lollipopPlot(maf_obj,gene="CCND3")

