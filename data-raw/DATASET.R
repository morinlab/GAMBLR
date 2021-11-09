## code to prepare `DATASET` dataset goes here

gene_blacklist = system.file("extdata", "gene_blacklist_with_IG.tsv", package = "GAMBLR") %>%
  read_tsv()

usethis::use_data(gene_blacklist, overwrite = TRUE)

grch37_all_gene_coordinates = system.file("extdata", "grch37_gene_coordinates.tsv", package = "GAMBLR") %>%
  read_tsv()

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

