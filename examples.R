package_path = "/projects/rmorin/projects/gambl-repos/gambl-rmorin/src/GAMBLR"
setwd(package_path)
devtools::install()
#set your working directory back to where you want to work
setwd("/projects/rmorin/projects/gambl-repos/gambl-rmorin")
library(GAMBLR)

# This code relies on a MySQL database at the GSC and will only work on that network, and only after you have configured your database connection (via ~/.my.cnf )
# See the Readme for more details


#get the metadata
my_meta = get_gambl_metadata()

#get all CN segments that overlap Chris' favourite gene

my_segments = get_cn_segments(chromosome="4",qstart=83274467,qend=83295149)

deleted_segments = my_segments %>% filter(log.ratio<0) #use some lower value if you want to be more stringent

annotated_segments = left_join(deleted_segments,my_meta,by=c("ID" = "sample_id"))

#Summarize by cohort

annotated_segments %>% pull(cohort) %>% table()

# CLL_GenomeCanada        DLBCL_ctDNA     DLBCL_Gascoyne DLBCL_GenomeCanada        DLBCL_HTMCP         DLBCL_ICGC
#         4                  2                  5                  9                  3                  5
# DLBCL_LSARP_Trios        DLBCL_Marra    FL_GenomeCanada            FL_ICGC          FL_Kridel    MALY_Other_ICGC
#       53                 10                  1                  4                  6                  3
# MM_mmsanger
#     4

#further restrict to just the first time point to avoid double-counting
annotated_segments %>% filter(time_point=="A") %>% pull(cohort) %>% table()

#get all mutations in the same gene
all_ssms = get_ssm_by_gene(gene_symbol=c("HNRNPD"),coding_only = TRUE)
all_ssms = all_ssms %>% as.data.frame()
#make a MAFtools object and plot a lollipop plot
library(maftools)
maf_obj = read.maf(all_ssms)
lollipopPlot(maf_obj,gene="HNRNPD")
