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

- [ ] I documented my function using [ROxygen style](https://jozef.io/r102-addin-roxytags/#:~:text=Inserting%20a%20skeleton%20%2D%20Do%20this,Shift%2BAlt%2BR%20).)

- [ ] All parameters for the function are described in the documentation and the function has a decriptive title. 

Example:
```
#' Use GISTIC2.0 scores output to reproduce maftools::chromoplot with more flexibility
#'
#' @param scores output file scores.gistic from the run of GISTIC2.0
#' @param genes_to_label optional. Provide a data frame of genes to label (if mutated). The first 3 columns must contain chromosome, start, and end coordinates. Another required column must contain gene names and be named `gene`. (truncated for example)
#' @param cutoff optional. Used to determine which regions to color as aberrant. Must be float in the range [0-1]. (truncated for example)
```

- [ ] My function uses a library that isn't already a dependency of GAMBLR and I made the package aware of this dependency using the function documentation `import` statment. 

Example:
```
#' @return nothing
#' @export
#' @import tidyverse ggrepel
```

## Checklist for changes to existing code

- [ ] I added/removed arguments to a function and updated documentation for all changed/new arguments

- [ ] I tested the new code for compatability with existing functionality in the Master branch (please provide a reprex of how you tested the original functionality)

