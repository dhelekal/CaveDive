# PLEASE SEE MANUSCRIPT FOR FULL DETAILS
## IN ORDER TO REPRODUCE RESULTS AND ANALYSES IN THE MANUSCRIPT FOLLOW THESE INTRUCTIONS

Note: Navigate to the directory containing this file
### First, run the batch simulations

1. MAKE SURE THE `run_mcmc` flag has been set to `TRUE` in the following files: `validate_zero_expansions.rmd`, `validate_one_expansion.rmd`, `validate_multiple_expansions.rmd`.
2. Render the above mentioned Rmarkdown files using `rmarkdown::render`. Note, this will take quite a while. 

### Second, run tree analysis notebooks for all GPSC phylogenies
This can be performed by executing `sh analyse_all_trees.sh Paper/GPSC_Trees` 

### Third, generate all other figures
1. MAKE SURE THE `run_mcmc` flag has been set to `TRUE` in the following files: `Paper/Figures/fig_uhlemann.R`, `Paper/Figures/figure2.R`.
2. Render all figures by running `sh make_figs.sh`

## Data references:
GPSC trees https://www.pneumogen.net/gps/
Uhlemann tree https://pathogen.watch/collection/pjycq3s0z83t-uhlemann-et-al-2014 NOTE: we only use the USA300 clade.
For full citations see manuscript.