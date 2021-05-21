# BrISCUT
**Br**eakpoint **I**dentification of **S**ignificant **C**ancer **U**ndiscovered **T**argets; detects genomic loci under selective pressures by analyzing partial SCNAs

Author: Juliann Shih, jshih@broadinstitute.org, juliann.shih@unlv.edu

Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu

BrISCUT (BrISCUT_pipeline.R) has 3 main steps:

1) Preprocessing step:
This takes a seg file and pumps out a folder named "breakpoint_files_[date]" with files similar to "PANCAN_3p_del_tel.txt" (found in "example" folder). The amplitude can be changed by changing "threshold".

2) Main BrISCUT script:
This takes the output of the preprocessing step and returns BrISCUT results.

3) Postprocessing and plotting:
BrISCUT will call combine_BrISCUT_results.py, plot_BrISCUT_results.R, plot_fig2.R, and calculate_effect_sizes.py which will aggregate significant peaks and plot within-cohort results as well as report peak effect sizes.

In the next few months we will be making this code base more user-friendly.

## How to run BrISCUT
### System requirements
This code has been tested on Macintosh and Linux using R version 3.6.1 and Python 2.7. The user must also install R packages ismev, extRemes, fitdistrplus, truncdist, segmented, parallel, dplyr, reticulate, ggplot2, pastecs, gridExtra, stringr, gtable, cowplot, and ggpubr.

### Running the script from R

```
source("BrISCUT_pipeline.R")
BrISCUT_pipeline(segfiles,ttlist,thedate,ci,infoloc,threshold)
  -segfiles [list of locations of segmented copy-number data; must correspond to ttlist]
  -ttlist [list of tumor types; must corresponds to segfiles]
  -thedate [enter date]
  -ci [confidence interval; suggested 0.95]
  -infoloc [location of file with start and end coordinates for each chromosome arm; example in docs]
  -threshold [log2 copy-number threshold for pSCNAs; suggested 0.2]
```


