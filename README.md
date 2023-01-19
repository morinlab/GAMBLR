# GAMBLR - an R package with convenience functions for working with GAMBL results

## Getting Started


This package relies heavily on a MariaDB/MySQL database hosted on the GSC network. Because of this, you can't use many of the functions outside of the network. You will also need a config file that contains the credentials that give you access to one of the databases. If you have the appropriate approvals, you can make a copy of the credentials file in the home directory of another GAMBL user.

```
cp /home/rmorin/.my.cnf ~/.my.cnf && chgrp morinlab ~/.my.cnf && chmod 750 ~/.my.cnf
```

It is extremely important that you ensure the file has the right unix group (`morinlab`) and is not readable by "all". The output of `ls -l` should look exactly as shown below (with your username shown instead of rmorin).

```
-rwxr-x--- 1 rmorin morinlab 54 Mar 31 08:41 /home/rmorin/.my.cnf
```

To confirm your connection works, connect to the database using the mysql client (installed on gphost03) using the following command:

```
mysql gambl_test
```

If you see a mysql prompt then you are ready to go. You can also see what tables are available and run queries here if you like.

```
show tables;
#don't forget the semicolon!
```

Instead of writing your own queries, it's more likely that you will want to use the dbplyr package and the convenience functions provided in GAMBLR. Adding new functionality is highly encouraged.

## Installation

Clone the repo to your home directory (not your gambl working directory)

```
git clone git@github.com:morinlab/GAMBLR.git
```

In Rstudio (on a gphost), set your working directory to the place you just cloned the repo.

```
setwd("~/GAMBLR")
```

Install the package in R by running the following command (requires the devtools package)

```
devtools::install()
```

I need to fix the dependencies to be auto-loaded. When you load the package it's a good plan to also load the dependencies (and install any that are missing).
```
library(GAMBLR)
library(dbplyr)
library(tidyverse)
library(data.table)
library(RMariaDB)
library(DBI)
library(ComplexHeatmap)
library(ggrepel)
library(ggthemes)
library(stats)
library(metaviz)
```

## Contributing

As GAMBL users (GAMBLRs, so to speak) rely on the functionality of this package, the Master branch is protected. All commits must be submitted via pull request on a branch. Please refer to the [GAMBL](https://github.com/morinlab/gambl#contribution-guidelines) documentation for details on how to do this.

# Run GAMBLR Remote On A Local Machine
This section is directed towards running GAMBLR remote on a *local machine*, for exploring a collection of GAMBLR functions. There are two different approaches to ge this to work. We will be going over both in the next sections.

## Approach 1 - Quick Start (SSH)
This section details how to deploy GAMBLR with limited functionalities, using ssh_session parameter. This approach requires either a working GSC VPN connection (or is directly accessible if connected to the GSC network).

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
3. Update the **paths** in your local [config.yml](https://github.com/morinlab/GAMBLR/blob/master/config.yml) to point to the recently cloned, local **gambl** folder (repo_base). In your favorite text editor, edit the line shown below (under *remote*).
```
repo_base: "path/to/your/local/gambl_repo/"
```
4. Set the **working directory** in Rstudio. Open Rstudio on your local machine and locate the repo you cloned in the step 2.
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
#### Run GAMBLR On Your Local Machine
9. Test if setup was successful (e.g call `get_gambl_metadata()` to retrieve meta data for all gambl samples).
```
get_gambl_metadata() %>%
  head()
```
## Approach 2 - The Full Instalation (Snakemake)
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
