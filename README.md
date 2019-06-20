# SamplingEffort
R code for the joint estimation of sampling effort and species distributions from species occurrences, a simulation example and application to Pl@ntNet data.

## Simulation experiment

### Prerequisites
R version 3.5.1 or superior (I didn't check earlier versions). 
Install librairies: glmnet, raster, ggplot2, data.table.

### Reproduce the experiment
0) Download all the repository (zip file __SamplingEffort-master.zip__) and unzip the whole directory __SamplingEffort-master__ locally.
1) Open script __Simu_and_graphs.R__ and modify the __dir__ variable to the location of the __SamplingEffort-master__ directory.
2) Run __Simu_and_graphs.R__ with R.
3) The individual graphs are saved as __.png__ images in the __SamplingEffort-master__ directory.

## Real data illustration

### Prerequisites
0) Having a machine with >60Gb of RAM. Required by glmnet during the fitting process (otherwise decreasing the number of background points per sampling cell - variable __n__ in __plantnet_effort.R__ - will approximately reduce the memory consumption by the same factor, but may entail higher bias in the model fit, or even identifiability problems).
1) Install librairies: glmnet, data.table, raster, rgdal, rgeos, grid, ggplot2.
2) Download all the repository (zip file __SamplingEffort-master.zip__) and unzip the whole directory __SamplingEffort-master__ locally.
3) Download the Pl@ntNet occurrences data:
Botella Christophe, Bonnet Pierre, Joly Alexis, Lombardo Jean-Christophe, & Affouard Antoine. (2019). Pl@ntNet queries 2017-2018 in France (Version 0) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.2634137
And move the files __PL_complete.csv__ and __taxaName_glc19SpId_InTest.csv__ to the __SamplingEffort-master__ directory.
4) Download the environmental rasters zip :
Botella Christophe. (2019). A compilation of environmental geographic rasters for SDM covering France (Version 1) [Data set]. Zenodo. http://doi.org/10.5281/zenodo.2635501
And unzip it inside the __SamplingEffort-master__ directory. 

### Reproduce the model fit and graphs 
1) Open script __plantnet_effort.R__ and modify the __dir__ variable to the location of the __SamplingEffort-master__ directory.
2) Run __plantnet_effort.R__ with R. (It will take a while... Go for day hike. If it crashes during glmnet fit, the )
3) 
