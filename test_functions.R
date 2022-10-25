# This is the working tests from some of the commonly used GAMBLR functions.
# Theoretically these should all work out-of-the-box and provide sensible outputs

library(GAMBLR)
library(tidyverse)
#starting point for most things: metadata for the samples you are working with (or all samples for one seq_type)
my_metadata = get_gambl_metadata() #by default, just the tumour genomes are returned

only_blgsp_metadata = get_gambl_metadata(case_set="BLGSP-study")
# set seq_type_filter to a different value if you need other seq_type (usually this will only be one at a time)
capture_metadata = get_gambl_metadata(seq_type_filter = 'capture')

non_duplicated_genome_and_capture = get_gambl_metadata(seq_type_filter=c('genome','capture'),seq_type_priority="genome")
#this gets you the metadata for both seq_types and drops sample_id that exist in both (picking genome as priority, in this case)


# best practices for getting SSMs
# if you want genome-wide calls for a handful of patients
patients = c("00-14595", "00-15201", "01-12047")
patients_maf = get_ssm_by_patients(these_patient_ids = patients, seq_type = "genome", subset_from_merge = FALSE)
patient_meta = get_gambl_metadata(seq_type_filter = "genome") %>% dplyr::filter(patient_id %in% patients)
patients_maf_2 = get_ssm_by_patients(these_samples_metadata = patient_meta,subset_from_merge = FALSE)

#or the same information but for a handful of samples
sample_ssms = get_ssm_by_samples(these_sample_ids=c("HTMCP-01-06-00485-01A-01D","14-35472_tumorA","14-35472_tumorB"))
hg38_ssms = get_ssm_by_samples(these_sample_ids=c("HTMCP-01-06-00485-01A-01D","14-35472_tumorA","14-35472_tumorB"),projection="hg38")

#You probably should never do this but this tests/shows its functionality
readr_sample_ssms = get_ssm_by_samples(these_sample_ids=c("HTMCP-01-06-00485-01A-01D","14-35472_tumorA","14-35472_tumorB"),subset_from_merge=TRUE,engine="readr")

slow_sample_ssms = get_ssm_by_samples(these_sample_ids=c("HTMCP-01-06-00485-01A-01D","14-35472_tumorA","14-35472_tumorB"),subset_from_merge=TRUE)


# you probably rarely need these functions but if you do want the same information for just one sample you can specify it with the sample_id or a single row from the metadata table containing your sample's information

this_sample_df = get_ssm_by_sample(this_sample_id = "HTMCP-01-06-00485-01A-01D", this_seq_type = "genome",tool_name = "slims-3", projection = "grch37")
capture_meta = get_gambl_metadata(seq_type_filter = "capture")
ssm_sample = get_ssm_by_sample(this_sample_id = "CASA0002_2015-03-10", projection = "grch37",augmented = T,these_samples_metadata = capture_meta)

#Coding mutations are probably where most people will want to start. This is how you get a MAF with just CDS variants (including Silent)
#There are some convenience options for requesting the result to be subset (filtered ) based on certain elements of the metadata
maf_data = get_coding_ssm(limit_cohort = c("BL_ICGC"), seq_type = "genome")
maf_data = get_coding_ssm(limit_samples = "HTMCP-01-06-00485-01A-01D", seq_type = "genome")

#To have ultimate control over what you get, I suggest you just subset the metadata the exact way you want then run it like so:
my_sample_metadata = get_gambl_metadata(seq_type_filter="capture") %>% dplyr::filter(pathology=="DLBCL")

maf_data = get_coding_ssm(these_samples_metadata = my_sample_metadata)

#mutations in regions of interest such as aSHM can be dealt with more efficiently using tabix-indexed MAFs
#WAY faster than loading the full file into RAM and filtering it

#for just one region:
my_mutations = get_ssm_by_region(region = "chr8:128,723,128-128,774,067")
#specifying chromosome, start and end individually instead of as a "region" string
my_mutations = get_ssm_by_region(chromosome = "8", qstart = 128723128, qend = 128774067)
#and in the familiar MAF format:
bcl2_all_details = get_ssm_by_region(region="chr18:60796500-60988073",basic_columns=T)

regions_bed = grch37_ashm_regions %>% mutate(name = paste(gene, region, sep = "_")) %>% head(20)
#this will get you the bare minimum basic 3-column information for mutated positions with non-MAF column naming
ashm_basic = get_ssm_by_regions(regions_bed = regions_bed)
#you can get all the normal MAF columns and header this way:
full_details_maf = get_ssm_by_regions(regions_bed = regions_bed, basic_columns=T)

#functions that use this to do fun stuff:
regions_bed = grch37_ashm_regions %>% mutate(name = paste(gene, region, sep = "_"))
matrix = get_ashm_count_matrix(regions_bed = head(regions_bed,15), seq_type="genome",these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(pathology=="DLBCL"))


#SVs
oncogenes_sv = get_combined_sv(oncogenes = c("MYC", "BCL2", "BCL6"))

#lazily get every SV in the table with default quality filters
all_sv = get_manta_sv()
#get all SVs for a single sample
some_sv = get_manta_sv(sample_id="94-15772_tumorA")
#get the SVs in a region around MYC
myc_locus_sv = get_manta_sv(region="8:128723128-128774067")


#Copy number
#Getting individual sets of segments
my_segments = get_cn_segments(region="chr8:128,723,128-128,774,067")
#specifying chromosome, start and end individually
my_segments = get_cn_segments(chromosome="8",qstart=128723128,qend=128774067)
#Asking for chromosome names to have a chr prefix (default is un-prefixed)
prefixed_segments = get_cn_segments(chromosome ="12",qstart = 122456912, qend = 122464036, with_chr_prefix = TRUE)

#this example uses the bundled bed file of regions containing lymphoma genes
#warning: This is pretty slow with the full bed file
cn_matrix = get_cn_states(regions_bed=head(grch37_lymphoma_genes_bed,15))

#this will just get the matrix for FL cases
cn_matrix = get_cn_states(these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(pathology=="FL"),all_cytobands = T,use_cytoband_name = T)
cn_matrix[cn_matrix>5]=5
