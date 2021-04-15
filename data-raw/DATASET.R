## code to prepare `DATASET` dataset goes here

grch37_oncogene = system.file("extdata", "oncogene_regions.grch37.tsv", package = "GAMBLR") %>%
  read_tsv(col_types="ciici")

usethis::use_data(grch37_oncogene, overwrite = TRUE)

grch37_partners = system.file("extdata", "superenhancer_regions.grch37.tsv", package = "GAMBLR") %>%
  read_tsv(col_types="ciici")

usethis::use_data(grch37_partners, overwrite = TRUE)

grch37_ashm_regions = system.file("extdata", "somatic_hypermutation_locations_GRCh37.txt", package = "GAMBLR") %>%
  read_tsv()

usethis::use_data(grch37_ashm_regions, overwrite = TRUE)

lymphoma_genes = system.file("extdata","lymphoma_genes.tsv",package="GAMBLR") %>%
  read_tsv(col_types="clllll")

usethis::use_data(lymphoma_genes, overwrite = TRUE)


