## code to prepare `DATASET` dataset goes here

grch37_oncogene = system.file("extdata", "oncogene_regions.grch37.tsv", package = "GAMBLR") %>%
  read_tsv(col_types="ciici")

usethis::use_data(grch37_oncogene, overwrite = TRUE)

grch37_partners = system.file("extdata", "superenhancer_regions.grch37.tsv", package = "GAMBLR") %>%
  read_tsv(col_types="ciici")

usethis::use_data(grch37_partners, overwrite = TRUE)

