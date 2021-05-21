# BrISCUT
**Br**eakpoint **I**dentification of **S**ignificant **C**ancer **U**ndiscovered **T**argets; detects genomic loci under selective pressures by analyzing partial SCNAs

Author: Juliann Shih, jshih@broadinstitute.org, juliann.shih@unlv.edu

Contact: Rameen Beroukhim, rameen_beroukhim@dfci.harvard.edu

BrISCUT has 3 main steps, encapsulated by these three scripts:

1) Preprocessing step (BrISCUT_preprocessing.py):
This takes a seg file (presently hard-coded into the do_arm function in lines 46 and 47) and pumps out a folder named "breakpoint_files_[date]" with files similar to "PANCAN_3p_del_tel.txt". The amplitude can be changed by changing "threshold". You'll have to play with a lot of the locations of files.

2) Main BrISCUT script (find_show_peaks_empirical_with_centromere_[date].R):
This takes the output of the preprocessing step and returns BrISCUT results.  main function (do_arm_gistic) 

3) Postprocessing and plotting:
BrISCUT will call combine_BrISCUT_results.py, plot_BrISCUT_results.R, and plot_fig2.R, which will aggregate significant peaks and plot within-cohort results.

## How to run BrISCUT
# System requirements
This code has been tested on Macintosh and Linux using R version 3.6.1 and Python 2.7.16. The user must also install R packages ismev, extRemes, fitdistrplus, truncdist, segmented, parallel, dplyr,  reticulate, ggplot2, pastecs, gridExtra, stringr, gtable, cowplot, and ggpubr.

# R



Present locations of the scripts in my folders:
/cga/meyerson/home/jshih/aneuploidy/BrISCUT_preprocessing.py
/cga/meyerson/home/jshih/aneuploidy/telomere/

Most recent example results are here: 
/cga/meyerson/home/jshih/aneuploidy/telomere/results_210319overlap_0.95
You can look in the summary folder under specific tumor types (e.g. "DLBC/summary/") for examples. The two files that you end up wanting are "[tt]_BrISCUT_results_cols_0.95.txt" and "[tt]_BrISCUT_results.txt", but there are also a few plots: "[tt]_BrISCUT_results_[neg/pos/negpos/amp/del].pdf"


