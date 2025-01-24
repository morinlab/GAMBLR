# <img src="img/figures/logo.png" align="right" alt="" width="120" />

![](https://github.com/morinlab/GAMBLR/actions/workflows/build_check.yml/badge.svg)

# GAMBLR - an R package with convenience functions for working with GAMBL results.

<details>

<summary>Which GAMBLR shoud I use?</summary>

If you are a member of the GAMBL consortium (Morin Lab or BC Cancer/CLC with approved access to restricted data), you came to the right place - please follow the instructions in the next section Installation.

If you are not a member of the GAMBL consortium (or waiting for your access to the restricted data) but want to take advantage of the powerful setup offered by GAMBL, please proceed to the non-restricted functionality of GAMBLR available through the [GAMBLR.open](https://github.com/morinlab/GAMBLR.open) repository. GAMBLR.open offers exact functionality of GAMBLR but with access to only published subset of data.

</details>

<details>

<summary>Installation</summary>

GAMBLR is an open-source package. It can be easily installed directly from GitHub:

```
devtools::install_github("morinlab/GAMBLR", repos = BiocManager::repositories())
```

This will install the full set of GAMBLR-verse children packages ([GAMBLR.data](https://github.com/morinlab/GAMBLR.data), [GAMBLR.helpers](https://github.com/morinlab/GAMBLR.results), [GAMBLR.utils](https://github.com/morinlab/GAMBLR.utils), [GAMBLR.viz](https://github.com/morinlab/GAMBLR.viz), [GAMBLR.results](https://github.com/morinlab/GAMBLR.results)) with all necessary dependencies. The latter child package ([GAMBLR.results](https://github.com/morinlab/GAMBLR.results)) requires access to the GSC resources and is not intended to be used outside of GSC. If you are interested in standalone functionality, please refer to the documentation of the [GAMBLR.open](https://github.com/morinlab/GAMBLR.open) package or any other individual child package.
</details>


<details>


<summary>Using GAMBLR</summary>

Once installed, the correct packages in correct order can be loaded with a regular `library` command:
```
library(GAMBLR)
```
This achieves sevaral goals at the same time: (i) the packages are loaded in correct order so there are no conflicts between functions existing in data and results packages; (ii) GAMBL-consortium users (Morin Lab members) have access to both public and restricted data; (iii) all children packages are loaded in one-liner instead of importing each child separately.
</details>


<details>

<summary>Contributing</summary>

If you have access to gphost, the easiest way to obtain and contribute to GAMBLR is to do this via cloning the repository

```
cd
git clone git@github.com:morinlab/GAMBLR.git
```

In your R editor of choice, set your working directory to the place you just cloned the repo.

```
setwd("~/GAMBLR")
```

Install the package in R by running the following command (requires the devtools package)

```
devtools::install()
```

As GAMBL users (GAMBLRs, so to speak) rely on the functionality of this package, the Master branch is protected. All commits must be submitted via pull request on a branch. Please refer to the [GAMBL](https://github.com/morinlab/gambl#contribution-guidelines) documentation for details on how to do this.

### New Functions

Please always ensure that the new function goes into the corresponding child package according to it's intended use. If you are not sure to which package the new function belongs to, please ask through opening new issue on this repository or starting new thread on Slack if you are the member of the Morin lab.

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

</details>
