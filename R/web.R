
#' @title Web Initialize GAMBL Site.
#'
#' @description Set up a fresh instance of a website to host on gitlab.
#'
#' @param site_base_name Base name for site.
#' @param base_directory Path to base directory.
#' @param my_name My name.
#' @param my_gitlab_email The email used for gitlab.
#'
#' @import workflowr
#' @export
#'
web_initialize_gambl_site = function(site_base_name,
                                     base_directory = "/home/rmorin/",
                                     my_name = "Ryan Morin",
                                     my_gitlab_email = "rdmorin@sfu.ca"){

  wflow_git_config(user.name = my_name, user.email = my_gitlab_email)
  setwd(base_directory)
  wflow_start(site_base_name)
  wflow_build()
}

web_add_update_page = function(path_to_markdown){
}
