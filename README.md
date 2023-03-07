![Build GAMBLR](https://github.com/morinlab/GAMBLR/actions/workflows/build_check.yml/badge.svg)

# GAMBLR - an R package with convenience functions for working with GAMBL results

## Installation
If you have access to gphost, the easiest way to obtain and run GAMBLR is to do this via Rstudio on a gphost. If you do not have access to gphost, please refer to the [Run Remote On A Local Machine](#Run Remote On A Local Machine) section. Assuming you are running Rstudio on gphost, clone the repo to your home directory (not your GAMBL working directory).

```
git clone git@github.com:morinlab/GAMBLR.git
```

In Rstudio (on a gphost), set your working directory to the place you just cloned the repo.

```
setwd("~/GAMBLR")
```

Install the package in R by running the following command (requires the devtools package)

```
devtools::install_github("morinlab/GAMBLR")
```

## Contributing

As GAMBL users (GAMBLRs, so to speak) rely on the functionality of this package, the Master branch is protected. All commits must be submitted via pull request on a branch. Please refer to the [GAMBL](https://github.com/morinlab/gambl#contribution-guidelines) documentation for details on how to do this.


#### Setup VPN Connection

1. You need a working GSC VPN connection to use this approach. For setting up a VPN connection see [this](https://www.bcgsc.ca/wiki/pages/viewpage.action?spaceKey=SysHelp&title=Learn+how+to+use+VPN) guide. Keep in mind that a **VPN connection is not needed** if your already connected to the GSC network.

#### Clone Repos, Update Paths, Install and Load R Packages


```
mkdir ~/git_repos
cd ~/git_repos #set as working directiory
git clone https://github.com/morinlab/gambl
git clone https://github.com/morinlab/GAMBLR
```

3. Update the **paths** in your local [config.yml](https://github.com/morinlab/GAMBLR/blob/master/config.yml) to point to the recently cloned, local **gambl** folder (repo_base). In your favorite text editor, edit the line shown below (under *remote*). Similarly, you will also need to edit the line above it to point to where you will eventually sync the GAMBL results.

```
remote:
    project_base: "/path/to/your/local/gambl_results_directory/"
    repo_base: "/path/to/your/local/gambl_repo/"
```

4. Set the **working directory** in Rstudio. Open Rstudio on your local machine and locate the repo you cloned previously.

```
setwd("~/git_repos/GAMBLR-master")
```

5. Install R packages to you local R studio session.

```
devtools::install()
install.packages("ssh")
```

6. Load packages.

```
library(GAMBLR)
library(ssh)
```

#### Setup SSH Connection in Rstudio
7. Specify your SSH connection in Rstudio. If your local machine username **does not** match your GSC username, you need to add *"your_username<!-- -->@gphost01.bcgsc.ca"* as the only argument to this function.

```
session = GAMBLR::get_ssh_session()
```

8. Execute the following in Rstudio console to make use of the updated paths in the [config.yml](https://github.com/morinlab/GAMBLR/blob/master/config.yml) from step 3.

```
Sys.setenv(R_CONFIG_ACTIVE = "remote")
```

#### Run GAMBLR
9. Test if setup was successful (e.g call `get_gambl_metadata()` to retrieve meta data for all gambl samples).

```
get_gambl_metadata() %>%
  head()
```

### Approach 2 - The Full Instalation (Snakemake)
This section details how to obtain GAMBLR with **full** functionality, using a dedicated snake file to retrieve all necessary files and dependencies.

#### Before You Get Started

1. Make sure you have a working SSH key **setup with a pass phrase**. If not, follow instructions at [GSC Wiki](https://www.bcgsc.ca/wiki/login.action?os_destination=%2Fpages%2Fviewpage.action%3FspaceKey%3DSysHelp%26title%3DSetup%2BPassword-less%2BSSH&permissionViolation=true). Warning, this will **not** work with a pass phrase-less SSH connection.

#### Clone Repos and Set Up Environment

2. Clone [GAMBL](https://github.com/morinlab/gambl) and [GAMBLR](https://github.com/morinlab/GAMBLR).

```
mkdir ~/git_repos
cd ~/git_repos
git clone https://github.com/morinlab/gambl
git clone https://github.com/morinlab/GAMBLR
```

3. On your local machine, make a new directory called **gambl_results**, for example.

```
mkdir ~/gambl_results/
```

4. Update paths under `remote` in your **local** [config.yml](https://github.com/morinlab/GAMBLR/blob/master/config.yml) (GAMBLR) to point to the recently cloned, local **gambl** folder (*repo_base*) and recently created **gambl_results** (*project_base*) folder. For example:

```
remote:
    project_base: "~/gambl_results/"
    repo_base: "~/git_repos/gambl-master/"
```

5. Copy the following files (from your recently cloned [GAMBLR](https://github.com/morinlab/GAMBLR) directory) into the folder from the previous step; `config.yml` and `get_gambl_results.smk`.

```
cp ~/git_repos/GAMBLR-master/config.yml ~/gambl_results/
cp ~/git_repos/GAMBLR-master/get_gambl_results.smk ~/gambl_results/
```

6. Add ENVVARS bash/zsh environment variables to your bashrc/zsh or some other way that will ensure they're in your session (e.g. you can set them manually each time if you want, just make sure they are set). For example in your local terminal run the following commands (with updated values...).

```
export GSC_USERNAME="your_gsc_username"
export GSC_KEY="path_to_SSH_key_with_passphrase_from_step_1"
export GSC_PASSPHRASE="passpharase_from_step_1"
```

#### Install GAMBLR In Local Rstudio

7. Open **Rstudio** (locally) and set the working directory to the folder you downloaded in step 2 (in the Rstudio console) and install GAMBLR.

```
setwd("~/git_repos/GAMBLR-master")
```

8. Install and load GAMBLR into your local R session.

```
devtools::install()
```

#### Create and Setup Snakemake Environment

9. In the terminal on your local machine, create a new snakemake environment from the [get_gambl_results.yml](https://github.com/morinlab/GAMBLR/blob/master/get_gambl_results.yml) file. Note that you can name this new environment whatever you would like. In this example, the new environment is called **snakemake_gambl**.

```
cd ~/gambl_results
conda env create --name snakemake_gambl --file ~/git_repos/GAMBLR-master/get_gambl_results.yml
```

10. Activate this newly created snakemake environment with:

```
conda activate snakemake_gambl
```

11. Retrieve necessary files (download a local copy of all files needed to run a collection of GAMBLR functions). It's strongly advised to use `--cores 1` for this, since it seems to be the more stable option. In addition, if your sync gets interrupted, you only need restart the syncing of 1 file, compared to if you run on multiple cores.

```
snakemake -s get_gambl_results.smk --cores 1
```

#### Use GAMBLR Functions Locally

12. In Rstudio (local), open [test_remote.R](https://github.com/morinlab/GAMBLR/blob/master/test_remote.R) in GAMBLR master folder.
13. Execute the following in Rstudio console to make use of the updated paths in the [config.yml](https://github.com/morinlab/GAMBLR/blob/master/config.yml) from step 5 (line 5 in [test_remote.r](https://github.com/morinlab/GAMBLR/blob/master/test_remote.R))

```
Sys.setenv(R_CONFIG_ACTIVE = "remote")
```

14. Check what files (if any) are currently missing.

```
check_gamblr_config()
```

15. You should now be all set to explore a collection of GAMBLR function remotely on your local machine. For example you could try the following test code to ensure your setup was successful. For a set of comprehensive examples and tutorials, please refer to the [test_remote.r](https://github.com/morinlab/GAMBLR/blob/master/test_remote.R) script.

```
get_gambl_metadata() %>%
  head()
```

**Note**, if your seeing the following message when trying to use GAMBLR, please ensure that the config/gambl repo is set up properly (step **5** and **13**) and/or remember to load the *remote* one (i.e `Sys.setenv(R_CONFIG_ACTIVE = "remote")`).

```
get_gambl_metadata(seq_type_filter = "capture") %>%
  pull(cohort) %>%
  table()

Error: '/projects/rmorin/projects/gambl-repos/gambl-rmorin/data/metadata/gambl_all_outcomes.tsv' does not exist.
```

## Contributing
As GAMBL users (GAMBLRs, so to speak) rely on the functionality of this package, the Master branch is protected. All commits must be submitted via pull request on a branch. Please refer to the [GAMBL](https://github.com/morinlab/gambl#contribution-guidelines) documentation for details on how to do this.

### For Developers
When designing new functions, please refer to guid-lines and best practices detailed [here](https://r-pkgs.org/). For your convenience, here is an empty function-skeleton that can be recycled when designing new GAMBLR functions. Ensure to always provide the required documentation for any new functions. See [this](https://r-pkgs.org/man.html#title-description-details) section for more details on best practices for documenting R functions. Unsure what information goes where in a function documentation? Here is a brief outline for what the different sections should include. For more information, see [this](https://r-pkgs.org/man.html#title).

#### Title
The title is taken from the first sentence. It should be written in sentence case, not end in a full stop, and be followed by a blank line. The title is shown in various function indexes (e.g. help(package = "somepackage")) and is what the user will usually see when browsing multiple functions.

#### Description
The description is taken from the next paragraph. It’s shown at the top of documentation and should briefly describe the most important features of the function.

#### Details
Additional details are anything after the description. Details are optional, but can be any length so are useful if you want to dig deep into some important aspect of the function. Note that, even though the details come right after the description in the introduction, they appear much later in rendered documentation.

```
#' @title
#'
#' @description
#'
#' @details
#'
#' @param a_parameter
#' @param another_parameter
#'
#' @return
#'
#' @import
#' @export
#'
#' @examples
#'

function_name = function(a_parameter,
                         another_parameter){
                         }
```
#### Example Function
For your convenience, as an example, here is a perfectly documented GAMBLR function, following the best practices detailed above.

```
#' @title ASHM Rainbow Plot
#'
#' @description Make a rainbow plot of all mutations in a region, ordered and coloured by metadata.
#'
#' @details This function creates a rainbow plot for all mutations in a region. Region can either be specified with the `region` parameter,
#' or the user can provide a maf that has already been subset to the region(s) of interest with `mutation_maf`.
#' As a third alternative, the regions can also be specified as a bed file with `bed`.
#' Lastly, this function has a variety of parameters that can be used to further customize the returned plot in many different ways.
#' Refer to the parameter descriptions, examples as well as the vignettes for more demonstrations how this function can be called.
#'
#' @param mutations_maf A data frame containing mutations (MAF format) within a region of interest (i.e. use the get_ssm_by_region).
#' @param metadata should be a data frame with sample_id as a column.
#' @param exclude_classifications Optional argument for excluding specific classifications from a metadeta file.
#' @param drop_unmutated Boolean argument for removing unmutated sample ids in mutated cases.
#' @param classification_column The name of the metadata column to use for ordering and colouring samples.
#' @param bed Optional data frame specifying the regions to annotate (required columns: start, end, name).
#' @param region Genomic region for plotting in bed format.
#' @param custom_colours Provide named vector (or named list of vectors) containing custom annotation colours if you do not want to use standartized pallette.
#' @param hide_ids Boolean argument, if TRUE, ids will be removed.
#'
#' @return ggplot2 object
#'
#' @import dplyr ggplot2
#' @export
#'
#' @examples
#' #basic usage
#' region = "chr6:90975034-91066134"
#' metadata = get_gambl_metadata()
#' plot = ashm_rainbow_plot(metadata = metadata, region = region)
#'
#' #advanced usages
#' mybed = data.frame(start=c(128806578,128805652,128748315), end=c(128806992,128809822,128748880), name=c("TSS","enhancer","MYC-e1"))
#' ashm_rainbow_plot(mutations_maf=my_mutations,metadata=my_metadata,bed=mybed)
#'
```
