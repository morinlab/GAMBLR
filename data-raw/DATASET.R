## code to prepare `DATASET` dataset goes here

gene_blacklist = system.file("extdata", "gene_blacklist_with_IG.tsv", package = "GAMBLR") %>%
  read_tsv()

usethis::use_data(gene_blacklist, overwrite = TRUE)

grch37_all_gene_coordinates = system.file("extdata", "grch37_gene_coordinates.tsv", package = "GAMBLR") %>%
  read_tsv() %>% dplyr::filter(grepl("PATCH",chromosome))

usethis::use_data(grch37_all_gene_coordinates, overwrite = TRUE)


hg38_oncogene = system.file("extdata", "oncogene_regions.hg38.tsv", package = "GAMBLR") %>%
  read_tsv(col_types="ciici")

usethis::use_data(hg38_oncogene, overwrite = TRUE)

grch37_oncogene = system.file("extdata", "oncogene_regions.grch37.tsv", package = "GAMBLR") %>%
  read_tsv(col_types="ciici")

usethis::use_data(grch37_oncogene, overwrite = TRUE)

hg38_partners = system.file("extdata","superenhancer_regions.hg38.tsv",package="GAMBLR") %>%
  read_tsv(col_types="ciici")

usethis::use_data(hg38_partners, overwrite = TRUE)

grch37_partners = system.file("extdata", "superenhancer_regions.grch37.tsv", package = "GAMBLR") %>%
  read_tsv(col_types="ciici")


grch37_ashm_regions = system.file("extdata", "somatic_hypermutation_locations_GRCh37.txt", package = "GAMBLR") %>%
  read_tsv() %>% mutate(name=paste0(gene,"-",region))

usethis::use_data(grch37_ashm_regions, overwrite = TRUE)

hg38_ashm_regions = system.file("extdata", "somatic_hypermutation_locations_GRCh38.txt", package = "GAMBLR") %>%
  read_tsv()

usethis::use_data(hg38_ashm_regions, overwrite = TRUE)

lymphoma_genes = system.file("extdata","lymphoma_genes.tsv",package="GAMBLR") %>%
  read_tsv(col_types="clllll")

usethis::use_data(lymphoma_genes, overwrite = TRUE)

#THis needs to be syncronized with the list above (it currently is out of sync!)
grch37_lymphoma_genes_bed = system.file("extdata","lymphoma_genes.grch37.bed",package="GAMBLR") %>%
  read_tsv()

usethis::use_data(grch37_lymphoma_genes_bed, overwrite = TRUE)

hg38_lymphoma_genes_bed = system.file("extdata","lymphoma_genes.hg38.bed",package="GAMBLR") %>%
  read_tsv()

usethis::use_data(hg38_lymphoma_genes_bed, overwrite = TRUE)


wright_genes_with_weights = system.file("extdata","WrightGenesWithWeights.txt",package="GAMBLR") %>%
  read.table(sep="\t",header=1) %>% rename(Ensembl_ID=EnsemblGeneID,Hugo_Symbol=GeneName)

usethis::use_data(wright_genes_with_weights, overwrite = TRUE)

dhitsig_genes_with_weights = system.file("extdata","DHITsigGenesWithWeights.txt",package="GAMBLR") %>%
  read.table(sep="\t",header=1) %>% rename(Ensembl_ID=ensembl_gene_id,Hugo_Symbol=GeneName)
usethis::use_data(dhitsig_genes_with_weights, overwrite = TRUE)

target_regions_hg38 = system.file("extdata","target_regions_hg38.txt",package="GAMBLR") %>%
  read.table(sep="\t",header=1)
usethis::use_data(target_regions_hg38, overwrite = TRUE)

target_regions_grch37 = system.file("extdata","target_regions_grch37.txt",package="GAMBLR") %>%
  read.table(sep="\t",header=1)
usethis::use_data(target_regions_grch37, overwrite = TRUE)

hotspot_regions_grch37 = system.file("extdata","hotspot_regions.grch37.tsv",package="GAMBLR") %>%
  read.table(sep="\t",header=1) %>%
  column_to_rownames("gene")
usethis::use_data(hotspot_regions_grch37, overwrite = TRUE)

hotspot_regions_hg38 = system.file("extdata","hotspot_regions.hg38.tsv",package="GAMBLR") %>%
  read.table(sep="\t",header=1) %>%
  column_to_rownames("gene")
usethis::use_data(hotspot_regions_hg38, overwrite = TRUE)

chromosome_arms_grch37 = system.file("extdata","chromosome_arms_grch37.tsv",package="GAMBLR") %>%
  read.table(sep="\t",header=1)
usethis::use_data(chromosome_arms_grch37, overwrite = TRUE)

chromosome_arms_hg38 = system.file("extdata","chromosome_arms_hg38.tsv",package="GAMBLR") %>%
  read.table(sep="\t",header=1)
usethis::use_data(chromosome_arms_hg38, overwrite = TRUE)


lymphgen_entrez = system.file("extdata","lymphgen_genes_entrez.txt",package="GAMBLR") %>%
  read_tsv()

entrez_map = system.file("extdata","hugo2entrez.tsv",package="GAMBLR") %>%
  read_tsv()

lymphgen_anno = left_join(lymphgen_entrez,entrez_map) %>% dplyr::rename("Hugo_Symbol"="Approved symbol") %>% dplyr::select(1:3)


library("biomaRt")
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl",version="hg38")
ensembl = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="grch37.ensembl.org", path="/biomart/martservice" ,dataset="hsapiens_gene_ensembl")
#need to get entrezgene_id, hgnc_symbol using ensembl_gene_id
gene_detail = getBM(attributes=c( 'ensembl_gene_id','entrezgene_id','hgnc_symbol'), 
      filters = 'ensembl_gene_id', 
      values = lymphoma_genes$ensembl_gene_id, 
      mart = ensembl,useCache = FALSE)

#fill in missing entrezgene_id

gene_detail[gene_detail$ensembl_gene_id =="ENSG00000036448","entrezgene_id"] = 9172

gene_detail[gene_detail$ensembl_gene_id =="ENSG00000065526","entrezgene_id"] = 23013

gene_detail[gene_detail$ensembl_gene_id =="ENSG00000102096","entrezgene_id"] = 11040

gene_detail[gene_detail$ensembl_gene_id =="ENSG00000172578","entrezgene_id"] = 89857

gene_detail[gene_detail$ensembl_gene_id =="ENSG00000205542","entrezgene_id"] = 7114

lymphoma_genes = left_join(lymphoma_genes,gene_detail) 
lymphoma_genes$LymphGen=FALSE
lymphoma_genes[lymphoma_genes$entrezgene_id %in% lymphgen_anno$NCBI_Gene_ID,"LymphGen"] = TRUE

#load the Reddy gene list
reddy_genes = system.file("extdata","reddy_genes.tsv",package="GAMBLR") %>%
  read_tsv() %>% dplyr::rename("hgnc_symbol"="Approved symbol")

usethis::use_data(reddy_genes, overwrite = TRUE)

lymphoma_genes$Reddy = FALSE
lymphoma_genes[lymphoma_genes$hgnc_symbol %in% reddy_genes$hgnc_symbol,"Reddy"]=TRUE

usethis::use_data(lymphoma_genes, overwrite = TRUE)

reddy_only = reddy_genes[which(!reddy_genes$hgnc_symbol %in% lymphoma_genes$hgnc_symbol),"hgnc_symbol"]

