# Pull Request Checklists

**Important:** When opening a pull request, keep only the applicable checklist and delete all other sections.

## Checklist for all PRs

### Required

- [ ] I tested the new code for my use case (please provide a reproducible example of how you tested the new functionality)

- [ ] I ensured all dplyr functions that commonly conflict with other packages are fully qualified. 

This can be checked and addressed by running `check_functions.pl` and responding to the prompts. Test your code _after_ you do this.

- [ ] I generated the documentation and checked for errors relating to the new function (e.g. `devtools::document()`) and added `NAMESPACE` and all other modified files in the root directory and under `man`. 

### Optional but preferred with PRs

- [ ] I updated and/or successfully knitted a vignette that relies on the modified code (which ones?)

## Checklist for New Functions

### Required

- [ ] I documented my function using [Roxygen style](https://jozef.io/r102-addin-roxytags/#:~:text=Inserting%20a%20skeleton%20%2D%20Do%20this,Shift%2BAlt%2BR%20).)

- [ ] Adequate function documentation (see [new-function documentation template](https://github.com/morinlab/GAMBLR#title) for more info)

- [ ] I have ran `devtools::document()` to add the newly created function to NAMESPACE (do not manually add anything to this file!).

Example:
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
#' @return ggplot2 object.
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
#' mybed = data.frame(start = c(128806578,
#'                              128805652,
#'                              128748315),
#'                    end = c(128806992,
#'                            128809822,
#'                            128748880),
#'                    name = c("TSS",
#'                             "enhancer",
#'                             "MYC-e1"))
#'
#' ashm_rainbow_plot(mutations_maf = my_mutations,
#'                   metadata = my_metadata,
#'                   bed = mybed)
#'
```

- [ ] My function uses a library that isn't already a dependency of GAMBLR and I made the package aware of this dependency using the function documentation `import` statement. 

Example:
```
#' @return nothing
#' @export
#' @import tidyverse ggrepel
```

## Checklist for changes to existing code

- [ ] I added/removed arguments to a function and updated documentation for all changed/new arguments

- [ ] I tested the new code for compatibility with existing functionality in the Master branch (please provide a reprex of how you tested the original functionality)

