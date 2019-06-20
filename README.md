# SamplingEffort
R code for the joint estimation of sampling effort and species distributions from species occurrences, a simulation example and application to Pl@ntNet data.

## Simulation experiment

### Prerequisites
R version 3.5.1 or superior (I didn't check earlier versions). 
Install librairies: glmnet, raster, ggplot2, data.table.

### Reproduce the experiment
0) Download the file __SamplingEffort-master.zip__ from the repository master branch and unzip it locally.
1) Open script __Simu_and_graphs.R__ and modify the __dir__ variable to the location of the directory __SamplingEffort-master__ (which contains the scripts just unzipped).
2) Run __Simu_and_graphs.R__ with R.

## Real data illustration

### Prerequisites
1) Install librairies: glmnet, raster, rgdal, rgeos, grid, ggplot2.
2) Download the Pl@ntNet occurrences data:
Botella Christophe, Bonnet Pierre, Joly Alexis, Lombardo Jean-Christophe, & Affouard Antoine. (2019). Pl@ntNet queries 2017-2018 in France (Version 0) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.2634137
And move the files PL_complete.csv and taxaName_glc19SpId_InTest.csv to the directory of script plantnet_effort.R.
3) Download the environmental rasters zip :
Botella Christophe. (2019). A compilation of environmental geographic rasters for SDM covering France (Version 1) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.2635501
And uncompress it to the directory of script plantnet_effort.R 
4) Make sure that the directory contains the file (should be contained

### Reproduce the model fit and graphs 
