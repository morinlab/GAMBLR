# This is the working tests from some of the commonly used GAMBLR functions.
# Theoretically these should all work out-of-the-box and provide sensible outputs

library(GAMBLR)
library(tidyverse)
#starting point for most things: metadata for the samples you are working with (or all samples for one seq_type)
my_metadata = get_gambl_metadata() #by default, just the tumour genomes are returned
only_blgsp_metadata = get_gambl_metadata(case_set = "BLGSP-study")

#set seq_type_filter to a different value if you need other seq_type (usually this will only be one at a time)
capture_metadata = get_gambl_metadata(seq_type_filter = "capture")

#this gets you the metadata for both seq_types and drops sample_id that exist in both (picking genome as priority, in this case)
non_duplicated_genome_and_capture = get_gambl_metadata(seq_type_filter = c('genome', 'capture'), seq_type_priority = "genome")


#utilities
#gene to region, for a single gene returned as "region"
myc_region = gene_to_region(gene_symbol = "MYC", genome_build = "grch37", return_as = "region")

#gene to region, for multiple genes (different genome_build) in a vector of characters, retuned as bed.
genes_bed = gene_to_region(gene_symbol = c("MYC", "BCL2"), genome_build = "hg38", return_as = "bed")

#gene to region, for multiple genes specified in an already subseted list of gene names, return format is data frame
fl_genes = dplyr::filter(lymphoma_genes, FL == TRUE) %>% pull(Gene)
fl_genes_df = gene_to_region(gene_symbol = fl_genes, genome_build = "grch37", return_as = "df")

#region to gene, simple example using a previously defined region (myc)
gene_in_this_region = region_to_gene(region = myc_region, gene_format = "hugo", genome_build = "grch37")

#region to gene, wider example using larger region (q-arm, chromosome 1) and using another geneom build (hg38) and gene format (ensembl)
chr1q_genes = region_to_gene(region = "chr1:142600000-249250621", gene_format = "ensembl", genome_build = "hg38")

#collate results, collate all genome samples, using cached results
genome_collated = collate_results(seq_type_filter = "genome", from_cache = TRUE)

#collate results, get collated results for all capture samples, using cached results
capture_collated = collate_results(seq_type_filter = "capture", from_cache = TRUE)

#collate results, use an already subset metadata table for getting collated results (cached)
fl_metadata = get_gambl_metadata(seq_type_filter = "genome") %>% dplyr::filter(pathology == "FL")
fl_collated = collate_results(seq_type_filter = "genome", join_with_full_metadata = TRUE, these_samples_metadata = fl_metadata, from_cache = TRUE)

#collate results, get collated results (cached) for all genome samples and join with full metadata
everything_collated = collate_results(seq_type_filter = "genome", from_cache = TRUE, join_with_full_metadata = TRUE)

#collate results, another example demonstrating correct usage of the sample_table parameter.
fl_samples = get_gambl_metadata(seq_type_filter = "genome") %>% dplyr::filter(pathology == "FL") %>% dplyr::select(sample_id, patient_id, biopsy_id)
fl_collated = collate_results(sample_table = fl_samples, seq_type_filter = "genome", from_cache = TRUE)

#regions to bin, bin chromosome 8 into 20000 bins
chr8q_bins = region_to_bins(chromosome = "8", start = 48100000, end = 146364022, bin_size = 20000)

#SSMs
#get ssm by patients, if you want genome-wide calls for a handful of patients
patients = c("00-14595", "00-15201", "01-12047")
patients_maf = get_ssm_by_patients(these_patient_ids = patients, seq_type = "genome", subset_from_merge = FALSE)

#get ssm by patients, you can also use an already subset metadata table directly (i.e these_samples_metadata)
patient_meta = get_gambl_metadata(seq_type_filter = "genome") %>% dplyr::filter(patient_id %in% patients)
patients_maf_2 = get_ssm_by_patients(these_samples_metadata = patient_meta,subset_from_merge = FALSE)

#get ssm by samples, or the same information but for a handful of samples
sample_ssms = get_ssm_by_samples(these_sample_ids = c("HTMCP-01-06-00485-01A-01D", "14-35472_tumorA", "14-35472_tumorB"))

#get ssm by samples, demonstrating another projection (hg38)
hg38_ssms = get_ssm_by_samples(these_sample_ids = c("HTMCP-01-06-00485-01A-01D", "14-35472_tumorA", "14-35472_tumorB"), projection = "hg38")

#get ssm by samples, you probably should never do this but this tests/shows its functionality
readr_sample_ssms = get_ssm_by_samples(these_sample_ids=c("HTMCP-01-06-00485-01A-01D","14-35472_tumorA","14-35472_tumorB"), subset_from_merge = TRUE, engine = "readr")
slow_sample_ssms = get_ssm_by_samples(these_sample_ids=c("HTMCP-01-06-00485-01A-01D","14-35472_tumorA","14-35472_tumorB"), subset_from_merge = TRUE)

#get ssm by sample, you probably rarely need these functions but if you do want the same information for just one sample you can specify it with the sample_id or a single row from the metadata table containing your sample's information
this_sample_df = get_ssm_by_sample(this_sample_id = "HTMCP-01-06-00485-01A-01D", this_seq_type = "genome", tool_name = "slims-3", projection = "grch37")

#get ssm by sample, demonstrating another seq type (capture)
capture_meta = get_gambl_metadata(seq_type_filter = "capture")
ssm_sample = get_ssm_by_sample(this_sample_id = "CASA0002_2015-03-10", projection = "grch37",augmented = TRUE, these_samples_metadata = capture_meta)

#get coding ssm, coding mutations are probably where most people will want to start. This is how you get a MAF with just CDS variants (including Silent)
#there are some convenience options for requesting the result to be subset (filtered ) based on certain elements of the metadata
maf_data = get_coding_ssm(limit_cohort = c("BL_ICGC"), seq_type = "genome")
maf_data = get_coding_ssm(limit_samples = "HTMCP-01-06-00485-01A-01D", seq_type = "genome")

#get coding ssm, to have ultimate control over what you get, I suggest you just subset the metadata the exact way you want then run it like so:
my_sample_metadata = get_gambl_metadata(seq_type_filter = "capture") %>% dplyr::filter(pathology == "DLBCL")
maf_data = get_coding_ssm(these_samples_metadata = my_sample_metadata)

#mutations in regions of interest such as aSHM can be dealt with more efficiently using tabix-indexed MAFs
#WAY faster than loading the full file into RAM and filtering it
#get ssm by region, for just one region:
my_mutations = get_ssm_by_region(region = "chr8:128,723,128-128,774,067")

#get ssm by region, specifying chromosome, start and end individually instead of as a "region" string
my_mutations = get_ssm_by_region(chromosome = "8", qstart = 128723128, qend = 128774067)

#get sm by region, and in the familiar MAF format:
bcl2_all_details = get_ssm_by_region(region = "chr18:60796500-60988073", basic_columns = TRUE)

#get ssm by region, this will get you the bare minimum basic 3-column information for mutated positions with non-MAF column naming
regions_bed = grch37_ashm_regions %>% mutate(name = paste(gene, region, sep = "_")) %>% head(20)
ashm_basic = get_ssm_by_regions(regions_bed = regions_bed)

#get ssm byr region, you can get all the normal MAF columns and header this way:
full_details_maf = get_ssm_by_regions(regions_bed = regions_bed, basic_columns = TRUE)

#get ashm count matrix, functions that use this to do fun stuff:
regions_bed = grch37_ashm_regions %>% mutate(name = paste(gene, region, sep = "_"))
matrix = get_ashm_count_matrix(regions_bed = head(regions_bed, 15), seq_type = "genome", these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(pathology == "DLBCL"))

#count ssm by region, example using minimal parameter as well as previously defined myc region
ssm_count_myc = count_ssm_by_region(region = myc_region, seq_type = "genome")

#count ssm by region, using chromosome, start and end parameters for defining the region of interest. Also subset to an already filtered metadata
ssm_counts = count_ssm_by_region(chromosome = 8, start = 128747680, end = 128753674, these_samples_metadata = fl_metadata, seq_type = "genome")

#count ssm by region, the last example for this function is making use of another parameter that allows the user to use an already called ssm object.
my_mutations = get_ssm_by_region(region = "chr8:128723128-128774067")
my_mutation_counts = count_ssm_by_region(all_mutations_in_these_regions = my_mutations, start = 128723128, end = 128774067, count_by = "sample_id", seq_type = "genome")

#get coding ssm status, tabulate mutation status for non-silent SSMs for a set of genes.
sample_ssms = get_ssm_by_samples(these_sample_ids=c("HTMCP-01-06-00485-01A-01D","14-35472_tumorA","14-35472_tumorB"))
coding_tabulated_df = get_coding_ssm_status(maf_data = sample_ssms, from_flatfile = TRUE, augmented = TRUE, seq_type = "genome", projection = "grch37", include_hotspots = FALSE, review_hotspots = FALSE, gene_symbols=c("MYC","KMT2D"))


#SVs
#get combined sv
oncogenes_sv = get_combined_sv(oncogenes = c("MYC", "BCL2", "BCL6"))

#get manta sv, lazily get every SV in the table with default quality filters
all_sv = get_manta_sv()

#get manta sv, get all SVs for a single sample
some_sv = get_manta_sv(sample_id = "94-15772_tumorA")

#get the SVs in a region around MYC
myc_locus_sv = get_manta_sv(region = "8:128723128-128774067")


#Copy number
#get cn segments, getting individual sets of segments
my_segments = get_cn_segments(region = "chr8:128,723,128-128,774,067", projection = "grch37", this_seq_type = "genome", streamlined = FALSE, from_flatfile = TRUE)

#get cn segments for specified region, using chromosome, qstart and qend parameters.
my_segments_2 = get_cn_segments(chromosome = "8", qstart = 128723128, qend = 128774067, projection = "hg38", this_seq_type = "genome", with_chr_prefix = TRUE, streamlined = FALSE, from_flatfile = TRUE)

#get cn segments, using a previously defined region (from gene_to_region) and filtering out all copy number states equal to 2
myc_cns = get_cn_segments(region = myc_region, projection = "grch37", this_seq_type = "genome", streamlined = FALSE, from_flatfile = TRUE, with_chr_prefix = FALSE) %>% dplyr::filter(CN != 2)

#get cn states, this example uses the bundled bed file of regions containing lymphoma genes, warning: This is pretty slow with the full bed file
cn_matrix = get_cn_states(regions_bed = head(grch37_lymphoma_genes_bed, 15))

#get cn states, this will just get the matrix for FL cases
cn_matrix = get_cn_states(these_samples_metadata = get_gambl_metadata() %>% dplyr::filter(pathology == "FL"), all_cytobands = TRUE, use_cytoband_name = TRUE)
cn_matrix[cn_matrix>5]=5

#get sample cn segments, return cn segments for one sample
sample_cn_segments = get_sample_cn_segments(this_sample_id = "HTMCP-01-06-00422-01A-01D", multiple_samples = FALSE, from_flatfile = TRUE, projection = "grch37", with_chr_prefix = FALSE, streamlined = FALSE)

#get sample cn segments, return cn segments for multiple samples
samples_cn_segments = get_sample_cn_segments(multiple_samples = TRUE, sample_list = c("00-15201_tumorA", "00-15201_tumorB"), from_flatfile = TRUE, projection = "hg38", with_chr_prefix = FALSE, streamlined = FALSE)

#assign cn to ssm, annotate mutations with copy number information
cn_list = assign_cn_to_ssm(this_sample = "HTMCP-01-06-00422-01A-01D", this_seq_type = "genome", coding_only = TRUE)

#get gene cn and expression, get copy number and expression for a single gene
MYC_cn_expression = get_gene_cn_and_expression(gene_symbol = "MYC")

#gene expression
#get gene expression, get gene expression for a single gene
metadata_mrna = get_gambl_metadata(seq_type_filter = "mrna", only_available = FALSE)
MYC_expr = get_gene_expression(hugo_symbols = c("MYC"), join_with = "mrna", metadata = metadata_mrna, all_genes = FALSE, from_flatfile = TRUE)

#get gene expression, another example, using ensembl IDs for multiple genes.
expression_data = get_gene_expression(ensembl_gene_ids = c("ENSG00000171791", "ENSG00000112149", "ENSG00000136997"), all_genes = FALSE, from_flatfile = TRUE)
