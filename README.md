# <img src="man/figures/logo.png" align="right" alt="" width="120" />
![](https://github.com/morinlab/GAMBLR/actions/workflows/build_check.yml/badge.svg)

# GABMLR - an R package with convenience functions for working with GAMBL results.
If you are viewing this page on Github, consider clicking [this link](https://morinlab.github.io/GAMBLR/) to go to the GAMBLR webpage and learn more about this package.

## Installation
If you have access to gphost, the easiest way to obtain and run GAMBLR is to do this via Rstudio on a gphost. If you do not have access to gphost, please refer to the [Run Remote On A Local Machine](#Run Remote On A Local Machine) section. Assuming you are running Rstudio on gphost, clone the repo to your home directory (not your GAMBL working directory).

```
git clone git@github.com:morinlab/GAMBLR.git
```

In Rstudio (on a gphost), set your working directory to the place you just cloned the repo.

```
setwd("~/GAMBLR-master")
```

Install the package in R by running the following command (requires the devtools package)

```
devtools::install()
```

## Running GAMBLR On Your Own Computer
If you don't have access to gphost on GSC, no worries, you can still execute GAMBLR functions in another way. Remote support was developed for this purpose. This section explains how to run GAMBLR remote on a *local machine* (i.e on your own computer). There are two different approaches to get this to work, both with its own advantages and limitations. We will be going over both in this next section.

### Approach 1 - Quick Start
This section details how to deploy GAMBLR with limited functionality. This approach requires either a working GSC VPN connection (or is directly accessible if connected to the GSC network).

#### Setup VPN Connection
1. You need a working GSC VPN connection to use this approach. For setting up a VPN connection see [this](https://www.bcgsc.ca/wiki/pages/viewpage.action?spaceKey=SysHelp&title=Learn+how+to+use+VPN) guide. Keep in mind that a **VPN connection is not needed** if your already connected to the GSC network.

#### Clone Repos, Update Paths, Install and Load R Packages
2. Clone [GAMBL](https://github.com/morinlab/gambl) and [GAMBLR](https://github.com/morinlab/GAMBLR) to your local computer. From your terminal run the following commands (folder structures can be whatever you want...)

```
mkdir ~/git_repos
cd ~/git_repos #set as working directiory
git clone https://github.com/morinlab/gambl
git clone https://github.com/morinlab/GAMBLR
```

3. Update the **paths** in your local [config.yml](https://github.com/morinlab/GAMBLR/blob/master/config.yml) (GAMBLR-master) to point to the recently cloned, local **gambl** folder (repo_base). In your favorite text editor, edit the line shown below (under *remote*). Similarly, you will also need to edit the line above it to point to where you will eventually sync the GAMBL results.

```
remote:
    project_base: "/path/to/your/local/gambl_results_directory/"
    repo_base: "/path/to/your/local/gambl_repo/"
```

4. Set the **working directory** in Rstudio. Open Rstudio on your local machine and locate the repo you cloned previously.

```
setwd("~/git_repos/GAMBLR-master")
```

5. Install GAMBLR in your local R studio.

```
devtools::install()
```

6. Load packages.

```
library(GAMBLR)
```

#### Set Config To Remote
7. Execute the following in Rstudio console to make use of the updated paths in the [config.yml](https://github.com/morinlab/GAMBLR/blob/master/config.yml) from step 3.

```
Sys.setenv(R_CONFIG_ACTIVE = "remote")
```

#### Run GAMBLR
8. Test if setup was successful (e.g call `get_gambl_metadata()` to retrieve meta data for all gambl samples).

```
get_gambl_metadata() %>%
  head()
```

### Approach 2 - The Full Installation (Snakemake)
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

4. Update paths under `remote` in your **local** [config.yml](https://github.com/morinlab/GAMBLR/blob/master/config.yml) (GAMBLR) to point to the recently cloned, local **gambl** folder (*repo_base*) and recently created **gambl_results** (*project_base*) folder. Also, update the *host* field to contain your username (you can use any gphost here). For example:

```
remote:
    project_base: "~/gambl_results/"
    repo_base: "~/git_repos/gambl-master/"
    ...
    host: "your_username@gphost01.bcgsc.ca"
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
9. In the terminal on your local machine, create a new snakemake environment from the [get_gambl_results.yml](https://github.com/morinlab/GAMBLR/blob/master/get_gambl_results.yml) file ([get_gambl_results_linux.yml](https://github.com/morinlab/GAMBLR/blob/master/get_gambl_results_linux.yml) for Linux). Note that you can name this new environment whatever you would like. In this example, the new environment is called **snakemake_gambl**.

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
12. In Rstudio (local), open [test_remote.R](https://github.com/morinlab/GAMBLR/blob/master/resources/test_remote.R) in GAMBLR master folder.
13. Execute the following in Rstudio console to make use of the updated paths in the [config.yml](https://github.com/morinlab/GAMBLR/blob/master/config.yml) from step 5 (line 5 in [test_remote.r](https://github.com/morinlab/GAMBLR/blob/master/resources/test_remote.R))

```
Sys.setenv(R_CONFIG_ACTIVE = "remote")
```

Alternatively, you can add the content of ~/git_repos/GAMBLR/.Rprofile to your ~/.Rprofile file. In this way, you do not need to enter the command above every time you start your R session.

```
cat ~/git_repos/GAMBLR/.Rprofile >> ~/.Rprofile
```

14. Check what files (if any) are currently missing.

```
check_gamblr_config()
```

15. You should now be all set to explore a collection of GAMBLR function remotely on your local machine. For example you could try the following test code to ensure your setup was successful. For a set of comprehensive examples and tutorials, please refer to the [test_remote.r](https://github.com/morinlab/GAMBLR/blob/master/resources/test_remote.R) script.

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

### New Functions
When designing new functions, please refer to guide-lines and best practices detailed [here](https://r-pkgs.org/). Ensure to always provide the required documentation for any new functions. See [this](https://r-pkgs.org/man.html#title-description-details) section for more details on best practices for documenting R functions. Unsure what information goes where in a function documentation? Here is a brief outline for what the different sections should include and as an example, [here](https://github.com/morinlab/GAMBLR/blob/master/R/viz.R#L2423) is an adequately documented GAMBLR function. For more information, see [this](https://r-pkgs.org/man.html#title).

#### Title
The title is taken from the first sentence. It should be written in sentence case, not end in a full stop, and be followed by a blank line. The title is shown in various function indexes (e.g. help(package = "some_package")) and is what the user will usually see when browsing multiple functions.

#### Description
The description is taken from the next paragraph. It’s shown at the top of documentation and should briefly describe the most important features of the function.

#### Details
Additional details are anything after the description. Details are optional, but can be any length so are useful if you want to dig deep into some important aspect of the function. Note that, even though the details come right after the description in the introduction, they appear much later in rendered documentation. If you want to add code in any other language other than R, this is also the sections to do so. For example, the new function relies on some bash code in order to utilize the GAMBLR code. You can detail such code here by simply adding a code block as you would in a regular markdown file.

#### Parameters
Detailed parameter descriptions should be included for all functions. Remember to state the required data types, default values, if the parameter is required or optional, etc.

#### Return
Specify the returned object, is it a data frame, a list, a vector or characters, etc.

#### Import
Always import all the packages from which you are calling any functions outside of base R and R [packages](https://cran.r-project.org/doc/FAQ/R-FAQ.html#R-Add_002dOn-Packages) that gets loaded per default. Remember to not import `tidyverse`, rather, import the individual packages from `tidyverse` that the function is depending on. If any packages that are not yet a part of the GAMBLR dependencies are needed for the function, the user needs to run `usethis::use_package("package_name")` in order to add any such new dependencies to DESCRIPTION. Warning, do not edit DESCRIPTION by hand, instead use the approach detailed here.

#### Export
Should this function be exported to NAMESPACE (i.e make it directly accessible for anyone who loads GAMBLR), or is the function considered to be an internal/helper function? In order to have the function populate NAMESPACE, the developer has to run `devtools::document()`. All functions that have the `@export` line in its documentation will be added to NAMESPACE. Helper functions should not include this in the function documentation. Note that such functions are still accessible with `GAMBLR:::helper_function_name`. If The new function is indeed a helper/internal function, ensure that this is made clear from both the function description and details (see [this](https://github.com/morinlab/GAMBLR/blob/master/R/utilities.R#L785) example). In addition, it should also be clear what purpose the helper function is serving (i.e what other GAMBLR functions are calling the helper function).

#### Examples
Please provide fully reproducible examples for the function. Ideally, the example should demonstrate basic usage, as well as more advanced usage with different parameter combinations. Note that examples can not extend over 100 characters per line, since this will cause the lines to  be truncated in the rendered PDF manual. In addition, the developer needs to load any packages (besides **GAMBLR**) that are needed to run the examples. For instance, if the example code calls `%>%`, `dplyr` or `magrittr` to make the pipe available for the example. It is advised to write your example in such a way that loading external packages are avoided as much as possible. Instead, prioritize base R as much as possible. In some cases, it is undesirable to have a function run its examples. This applies to functions that are writing files and helper functions. To avoid any such examples to run, simply wrap the example in:
```
\dontrun{
do_not_run = some_function()
}
```

#### Helper Function Specific Instructions
If the newly added function is a utilized by GAMBLR as an internal or helper function, you should also add the `@noRd` field to the function documentation. This prevents the function to have an `.Rd` file created and populated in the `man/` folder. This is important, since such functions should not be represented on the website that is being built from the source code. In addition, make sure that you also followed the helper function specific instructions under **Examples** and **Export**.

### Testing New Functions
So you have added a new function (carefully following the steps in the previous section!) and you are obviously extremely proud and eager to test it out (and let others test it). There are basically two different approaches to do so. 

#### Option 1
Your first option, and likely the preferred route to take, is to make sure that the working directory in R studio is set to the GAMBLR folder with your updated code and then run `devtools::load_all()` to load all the functions available in the `R/` folder of thee same repo. This should make all such functions available to call.

#### Option 2
As an alternative, you can also run `devtools::install()` from the updated GAMBLR directory. As the name implies, this will install the complete package complete with dependencies, remotes, etc. **Note**, if you run with the second option, make sure to restart your R session with `.rs.restartR()` after installing the package and then load GAMBLR with `library(GAMBLR)`. Now you have installed the updated branch of GAMBLR and are free to call any functions available in the `R/` 

### Function Documentation Template
For your convenience, here is a documentation template for GAMBLR functions.
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
#' #this is an example
#' ###For your reference, this line is exactly 100 characters. Do not exceed 100 characters per line
#'
function_name = function(a_parameter,
                         another_parameter){
                         }
```
