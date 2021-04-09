## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup,warning=FALSE,message=FALSE----------------------------------------
# Load the GAMBLR package and other packages required by these examples. 
library(GAMBLR)
require(tidyverse)
require(maftools)
require(circlize)

## ----connect_show_tables,message=FALSE,warning=FALSE--------------------------
con <- DBI::dbConnect(RMariaDB::MariaDB(), dbname = "gambl_test")
# use DBI function to list the tables
all_table_names = DBI::dbListTables(con)
print(all_table_names)
# peek at the contents of each table. 
for(table_name in all_table_names){
  table_db <- tbl(con, table_name)
  print(table_name)
  
}


## ----metadata,message=FALSE,warning=FALSE-------------------------------------

my_metadata = get_gambl_metadata()

# reduce to some of the more useful metadata fields

my_metadata = my_metadata %>% select(sample_id,biopsy_id,myc_ba,cohort,pathology)

print(head(my_metadata))


## ----mutation_totals, out.width="75%"-----------------------------------------
hg19_maf = tbl(con,"maf_slms3_hg19_icgc")
mutation_counts = hg19_maf %>% group_by(Tumor_Sample_Barcode) %>% tally()
print(head(mutation_counts))

mutation_counts  %>% ggplot() + geom_histogram(aes(x=n),bins=100) + xlim(c(10,25000))


## ----join_metadata------------------------------------------------------------
sample_meta = tbl(con,"sample_metadata") %>% filter(seq_type == "genome" & tissue_status == "tumour")
#if we only care about genomes, we can drop/filter anything that isn't a tumour genome
#The key for joining this table to the mutation information is to use sample_id. Think of this as equivalent to a library_id. It will differ depending on what assay was done to the sample. 

biopsy_meta = tbl(con,"biopsy_metadata") %>% select(-patient_id) %>% select(-pathology) %>% select(-time_point) %>% select(-EBV_status_inf)
# this table is keyed on biopsy_id. One biopsy can have more than one sample_id. This table must be joined to the sample table using biopsy_id. Because of some redundancy in columns between the two tables, I've dropped all redundant columns with the exception of biopsy_id.

all_meta = left_join(sample_meta,biopsy_meta,by="biopsy_id")
# IMPORTANT: the dbplyr package uses mysql queries under the hood to lazily retrieve the datat you need for each table on-the-fly. To properly use the efficiency of indexing and joins in MySQL, don't convert your tables into data frames until you're done all the necessary joins. 

# for working with the metadata when querying mutations, you'll need another join. 
# left join the metadata to the mutations. Somewhat confusingly, the join is on columns with different names. To deal with this, use: by = c("Tumor_Sample_Barcode" = "sample_id")



## ----manta,out.width="80%"----------------------------------------------------

myc_locus_sv = get_manta_sv(region="8:128723128-128774067",pass=FALSE,with_chr_prefix=TRUE) 

# we can override default that requires SV to meet the Manta "Pass" filtering criterion
# here we are also asking for the chromosomes to be named with a chr prefix (for circlize compatability)

bed1 = myc_locus_sv %>% select(CHROM_A,START_A,END_A,tumour_sample_id)
bed2 = myc_locus_sv %>% select(CHROM_B,START_B,END_B,tumour_sample_id)
colnames(bed1)=c("chrom","start","end","sample_id")
colnames(bed2)=c("chrom","start","end","sample_id")

myc_locus_sv = myc_locus_sv %>% select(tumour_sample_id,CHROM_A,START_A,CHROM_B,START_B,STRAND_A,STRAND_B,VAF_tumour)

print(head(myc_locus_sv))

circos.initializeWithIdeogram()
circos.genomicLink(bed1, bed2,col = rand_color(nrow(bed1),transparency=0.5))


## ----rainfall_plots,out.width="100%", fig.dim=c(8,3)--------------------------
# we can also query the database to get a MAF per patient for patient-centric visualizations. 
# note that here I'm not restricting to only coding variants
example_dlbcl = hg19_maf %>% filter(Tumor_Sample_Barcode == "HTMCP-01-06-00422-01A-01D")
example_dlbcl_df = example_dlbcl %>% as.data.frame()
example_dlbcl_maf = read.maf(example_dlbcl_df)
rainfallPlot(example_dlbcl_maf)

## ----rainbow_plot_ashm, out.width="100%", fig.dim=c(8,3)----------------------

# set up some coordinates to annotate in your plot (optional)
mybed = data.frame(start=c(128806578,128805652,128748315), end=c(128806992,128809822,128748880), name=c("TSS","enhancer","MYC-e1"))

# get the mutations within a region of interest
my_mutations = get_ssm_by_region(region="chr8:128,743,606-128,820,015")


ashm_rainbow_plot(mutations_maf=my_mutations,metadata=my_metadata,bed=mybed)

## ----get_seg_data-------------------------------------------------------------

my_segments = get_cn_segments(chromosome="4",qstart=83274467,qend=83295149)
print(head(my_segments))

deleted_segments = my_segments %>% filter(log.ratio<0) #use some lower value if you want to be more stringent

annotated_segments = left_join(deleted_segments,my_metadata,by=c("ID" = "sample_id"))

annotated_segments %>% pull(cohort) %>% table()



## ----ssms_and_maftools,out.width="100%", fig.dim=c(8,3)-----------------------

all_ssms = get_ssm_by_gene(gene_symbol=c("CCND3"),coding_only = TRUE)

all_ssms = all_ssms %>% as.data.frame()

# make a MAFtools object and plot a lollipop plot

maf_obj = read.maf(all_ssms)
lollipopPlot(maf_obj,gene="CCND3")

