# InfantFeedingAnalysis
Analysis of infant feeding practices by maternal HIV status reported in Demographic and Health Surveys in sub-Saharan Africa

# Description
The Demographic and Health Survey (DHS) Program collects population survey data from over 90 countries. In this repository, we use data from DHS Program surveys in sub-Saharan Africa to estimate the duration of breastfeeding by region (Central, Eastern, Southern, or Western Africa) by HIV status. The output of this script is intended for use in the [Spectrum software suite](https://avenirhealth.org/software-spectrum.php) many countries use to estimate HIV burden.

# Requirements
This analysis is implemented in R and uses the [rdhs](https://cran.r-project.org/web/packages/rdhs/index.html) package to access survey microdata from the [DHS Program](https://dhsprogram.com/). You must register an account with the DHS Program and request access to survey microdata, including HIV testing datasets, to use the code in this repository. Once you are authorized to access survey datasets, you will need to configure rdhs to access survey data (see the [rdhs documentation](https://github.com/ropensci/rdhs#readme) for details).

# How to use this code
The file bf-model.R is the main entry point for this code. You can run this from R by pointing your current working directory at this repository, then running `source(bf-model.R)`. This will take some time, as the program will first download survey datasets from the DHSProgram, perform some data reorganization, then use [rstan](https://cran.r-project.org/web/packages/rstan/index.html) to fit a statistical model to survey data. Once the model has been fitted, the script will produce figures showing the model fit to data from each survey. After fitting is completed, you can run `source(gen-figure.R)` to plot regional average estimates for 2005, 2010, and 2015.
