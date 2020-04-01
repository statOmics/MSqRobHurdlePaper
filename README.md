# MSqRob leaps over the hurdle: uniting peptide counts and peak-intensity-based quantitative proteomics

This repository contains all required code to reproduce the analyses for our pulication which has just been accepted in Analytical Chemistry (see https://pubs.acs.org/doi/pdf/10.1021/acs.analchem.9b04375). 

## Overview of the repository

This repository contains the following main folders:

- R: this folder stores the main functions to execute MSqRob (functions.R), the quasibinomial model (functions_Hurdle.R), the hurdle model (functions.R and functions_Hurdle.R), some plotting functions (PlotFDRTPR.R) and a bug fix for MSstats so that it is able to import MaxQuant data (MaxQtoMSstatsFormat.R). If you want to know how are methods are implemented exactly, this is the place to be.

- datasets: this folder contains all search result files for the CPTAC dataset and the HEART dataset and some extra data which we added ourselves, such as missing gene names in the HEART dataset.

- analyses: this folder contains all the files that were used to (1) analyse the CPTAC dataset, (2) analyse the HEART dataset and (3) make the plots in the main article and the plots and tables in supplementary material.

- save_files_missing_PRIDE, save_files_CPTAC and save_files_PXD006675 contain saved intermediary and final results. These come in very handy to quickly reproduce e.g. a plot without (re-)running all slow parts of the code.

## How to reproduce our data

All analyses were done in RMarkdown (https://rmarkdown.rstudio.com). If you want to reproduce the html files, you will need an installation of R (https://www.r-project.org) and RStudio (https://www.rstudio.com/products/rstudio/download/) or the packages needed to work with RMarkdown. Next, download this project to your computer. If you don't want to reproduce the html files, you can also copy the R code chuncks in the .Rmd files to R and execute them.

### Analysis of the CPTAC data

Run the files in analyses/CPTAC. Many temporary results have been saved and are loaded by default to speed up the workflow. However, you can always check the code that was used to generate them. Thanks to the saved files, each file is created in such a way that it can be run independently of the other files.

Files are numbered so that if you don't want to make use of the saved intermediary results, you can run them in ascending order just like we did to generate the results.

### Analysis of the HEART data

**Important:** for the PXD006675 (HEART) dataset, you will first need to download the peptides.txt and proteinGroups.txt file from the PRIDE repository at https://www.ebi.ac.uk/pride/archive/projects/PXD006675/files. Download the file "search.zip" and extract it. The peptides.txt and proteinGroups.txt files will be in the "txt" folder. Place these files under "datasets/PXD006675".

Another file that could not be made available because of it size is "res_HEART.RData". Hence, the file "1_analysis_PXD006675_MSqRob.Rmd" won't run by default. Under "2.3. Fit robust mixed models and test contrasts", you will need to uncomment and execute the following piece of code:

`res.HEART <- do_mm(formula = form, msnset = peptides.HEART3, group_var = gene.name, contrasts = contrasts, type_df = "conservative", max_iter = 20) #Lasts about 1 hour on our system`

And then save the result to your system as follows:

`save(res.HEART, file = "res_HEART.RData")`

Then place the res_HEART.RData file under the "save_files_PXD006675" directory. Then, "1_analysis_PXD006675_MSqRob.Rmd" will run as intended. Alternatively, if you don't want to execute MSqRob, you can skip this part and immediately load the res.HEART.full object (last line of code under "2.7 Remove contrasts with only one sample identified per condition").

Each file is again numbered so that you can run them in ascending order if you don't want to make use of the saved intermediary results.

### Making the plots

Again, thanks to the saved files, each of these files will rund independently. "1_comparison_plots.Rmd" creates all the plots in the main article, while "2_supp_comparison_plots.Rmd" creates all the plots in supplementary information. "2_supp_comparison_plots.Rmd" also outputs the data on which the supplementary tables are based.

## Reproduce our bioRxiv preprint

All analyses for our bioRxiv preprint at https://www.biorxiv.org/content/10.1101/782466v1 can be found at the corresponding GitHub branch: https://github.com/statOmics/MSqRobHurdlePaper/tree/bioRxiv.

## Help

If anything is unclear or doesn't work, please do not hesitate to contact ludger.goeminne@epfl.ch or raise an issue on GitHub.

## Citation

If you make use of the functions in this repository, please refer to: Ludger J.E. Goeminne, Adriaan Sticker, Lennart Martens, Kris Gevaert and Lieven Clement (2020). **MSqRob takes the missing hurdle: uniting intensity- and count-based proteomics**. *Analytical Chemistry*. doi: https://doi.org/10.1021/acs.analchem.9b04375.
