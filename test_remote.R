
setwd("~/git/gambl/") # set this path to point to your local clone of the gambl repo
# in your config.yml file set repo_base under the remote config (at the top) to point to this path also
library(GAMBLR)
Sys.setenv(R_CONFIG_ACTIVE= "remote")


config::get("project_base")
#[1] "/Users/rmorin/gambl_results/"

all_meta = get_gambl_metadata()

collated = collate_results(from_cache=TRUE)

collated = collate_results(from_cache=TRUE,seq_type_filter="capture")

