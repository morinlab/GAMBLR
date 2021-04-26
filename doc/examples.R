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
require(data.table)
require(rtracklayer)
require(RMariaDB)
require(DBI)

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

#An important feature for reproducibility is that we all use the same exact subset of samples when performing various level-3 analyses for a study. This function effectively "locks in" a set of cases for a study based on some filters it applies automatically based on a study set identifier. For retrieving the metadata for all BLGSP cases with WGS data, you can use the following:

blgsp_metadata = get_gambl_metadata(case_set="BLGSP-study")

#as you can see this spans several cohorts and doesn't even just include BL pathology (for complex reasons)

blgsp_metadata %>% pull(cohort) %>% table()

#If you want the sample_id (i.e. the Tumor_Sample_Barcode) for all these cases, for example to subset a MAF file, you can extract them into a vector:

blgsp_sample_ids = pull(blgsp_metadata,sample_id)
length(blgsp_sample_ids)


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


## ----oncogene_annotations, out.width="100%"-----------------------------------

unannotated_sv = get_manta_sv() 

#in this example, let's just look at the SVs annotated as likely driving BCL6 expression
annotated_sv = annotate_sv(unannotated_sv,with_chr_prefix = TRUE) %>% 
    filter(!is.na(partner)) %>% 
    #select(-entrez) %>% 
  filter(gene=="BCL6") %>% 
  as.data.frame()

print(head(annotated_sv))

#for labelling, get the unique set of partners
to_label =unique(annotated_sv$partner)
partner_label = grch37_partners %>% 
  filter(gene %in% to_label) %>% mutate(chrom = paste0("chr",chrom)) %>%
  as.data.frame()

bed1 = annotated_sv %>% select(chrom1,start1,end1,tumour_sample_id,fusion)
bed2 = annotated_sv %>% select(chrom2,start2,end2,tumour_sample_id,fusion)
colnames(bed1)=c("chrom","start","end","sample_id","fusion")
colnames(bed2)=c("chrom","start","end","sample_id","fusion")

#circos.initializeWithIdeogram()
circos.clear()
circos.initializeWithIdeogram(plotType = NULL,chromosome.index = paste0("chr", c(1:22,"X")))

circos.genomicLabels(partner_label, labels.column = 4, side = "outside", cex=0.4,col='black')

circos.genomicIdeogram()
circos.genomicLink(bed1, bed2,col = rand_color(nrow(bed1),transparency=0.5))



## ----oncogene_annotations_myc, out.width="100%"-------------------------------

unannotated_sv = get_manta_sv() 

#Now let's try the SVs annotated as likely driving MYC expression
annotated_sv = annotate_sv(unannotated_sv,with_chr_prefix = TRUE) %>% 
    filter(!is.na(partner)) %>% 
  filter(gene=="MYC") %>% 
  as.data.frame()

print(head(annotated_sv))

#for labelling, get the unique set of partners
to_label =unique(annotated_sv$partner)
partner_label = grch37_partners %>% 
  filter(gene %in% to_label) %>% mutate(chrom = paste0("chr",chrom)) %>%
  as.data.frame()
onco_label = grch37_oncogene %>% mutate(chrom = paste0("chr",chrom)) %>%
  as.data.frame()

bed1 = annotated_sv %>% select(chrom1,start1,end1,tumour_sample_id,fusion)
bed2 = annotated_sv %>% select(chrom2,start2,end2,tumour_sample_id,fusion)
colnames(bed1)=c("chrom","start","end","sample_id","fusion")
colnames(bed2)=c("chrom","start","end","sample_id","fusion")

circos.initializeWithIdeogram(plotType = NULL,chromosome.index = paste0("chr", c(1:22,"X")))

#circos.genomicLabels(partner_label, labels.column = 4, side = "outside", cex=0.4,col='black')

#circos.genomicLabels(onco_label, labels.column = 4, side = "outside", cex=0.4,col='red')

all_labels = rbind(onco_label,partner_label)
cols= c(rep("red",length(onco_label$chrom)),rep("black",length(partner_label$chrom)))
circos.genomicLabels(all_labels, labels.column = 4, side = "outside", cex=0.7,col=cols)
circos.genomicIdeogram()
circos.genomicLink(bed1, bed2,col = rand_color(nrow(bed1),transparency=0.5))



## ----lifting_over-------------------------------------------------------------
bedpe_hg19 = unannotated_sv %>% head(20) %>% as.data.frame()
print(head(bedpe_hg19))
bedpe_hg38 = liftover_bedpe(bedpe_df = bedpe_hg19,target_build = "hg38")
print(head(bedpe_hg38))


## ----rainfall_plots,out.width="100%", fig.dim=c(8,3)--------------------------
# we can also directly query the database to get a MAF per patient for patient-centric visualizations. # this is not using the GAMBLR functions but shows an example of how unsupported queries can be accomplished. Beware queries that will return many thousands of variants. Thes will be slow and may fail if they're too greedy 
# note that here I'm not restricting to only coding variants
example_dlbcl = hg19_maf %>% filter(Tumor_Sample_Barcode == "HTMCP-01-06-00422-01A-01D")
example_dlbcl_df = example_dlbcl %>% as.data.frame()
example_dlbcl_maf = read.maf(example_dlbcl_df)
rainfallPlot(example_dlbcl_maf)
DBI::dbDisconnect(con)

## ----copy_number_and_vaf------------------------------------------------------
#use the same sample as the previous example
my_sample = "HTMCP-01-06-00422-01A-01D"

copy_number_vaf_plot(this_sample=my_sample)

#what if we want to focus on putative driver mutations? You can restrict this plot just to coding mutations and label genes of your choice.

#use the built in lymphoma gene list and subset for BL or DLBCL

my_genes=lymphoma_genes %>% filter(BL==TRUE | DLBCL == TRUE) %>% pull(Gene)

copy_number_vaf_plot(this_sample=my_sample,coding_only = TRUE,genes_to_label = my_genes)




## ----rainbow_plot_ashm, out.width="100%", fig.dim=c(8,3)----------------------

# set up some coordinates to annotate in your plot (optional)
mybed = data.frame(start=c(128806578,128805652,128748315), end=c(128806992,128809822,128748880), name=c("TSS","enhancer","MYC-e1"))

# get the mutations within a region of interest
# note that we can specify the query chromosome with or without a chr prefix and it will be handled elegantly
my_mutations = get_ssm_by_region(chromosome="chr8",qstart=128743606,qend=128820015)

ashm_rainbow_plot(mutations_maf=my_mutations,metadata=my_metadata,bed=mybed)

## ----more_ashm_plotting-------------------------------------------------------
# Handy function that provides a vector of colours for giving points for different pathology/subgroups reproducible and distinguishable colours
lymphgen_colours = get_gambl_colours(classification="lymphgen")
# This package comes with some custom (curated) data such as the regions recurrently affected by hypermutation in B-NHLs
ashm_multi_rainbow_plot(regions_to_display=c("BCL2-TSS","MYC-TSS","SGK1-TSS","IGL"),custom_colours = lymphgen_colours)


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

## ----get_big_maf,width="100%", fig.dim=c(8,5)---------------------------------
# load the master merged MAF (coding only). It's usually more efficient to do this than to try to add filters to this query. Instead, just filter the data afterward
maf_data = get_coding_ssm() 
maf_metadata = get_gambl_metadata() %>% dplyr::rename("Tumor_Sample_Barcode"="sample_id")
maf = read.maf(maf_data,clinicalData = maf_metadata)
# subset metadata, for example, to only BLGSP samples
blgsp_metadata = get_gambl_metadata(case_set="BLGSP-study")
# subset maf
blgsp_maf <- subsetMaf(maf, tsb=blgsp_sample_ids)
# perform analysis for subset maf
oncoplot(blgsp_maf,clinicalFeatures = c("sex","cohort"),sortByAnnotation = TRUE,genesToIgnore = c("TTN","LILRB1"))

