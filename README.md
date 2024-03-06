# Core repository - Assessment of the detrimental effects of collinearity in classical and transformation models #

This repository contains the core files of the thesis of Jerome Sepin for the Master Program in Biostatistics at the University of Zurich.
The repository with the files to execute the simulation study and all the results can be found in https://bitbucket.org/jsepin/simulation/src/master .

### How do I work with these repositories? ###
* Clone the following two git repositories into the same directory: https://bitbucket.org/jsepin/STA495MT_JS/src/master (**STA495MT_JS**) and https://bitbucket.org/jsepin/simulation/src/master (**simulation**) or download them manually.
* Compile (*Build All*) the **STA495MT_JS/STA495MasterThesis/report/report.Rproj** project. This will provide the parameters for the experimental conditions used in the simulation.
* Run the **simulation/simulation_total.R** file. This will perform the whole simulation and takes about 12h. You need the experimental conditions which get saved in **STA495MT_JS/STA495MasterThesis/data/boston_parameters.rds** .
* Compile (*Build All*) the simulation/results_simulation/results_simulation.Rproj project. This will provide the figures for the results and the demonstration data frame (**simulation/data/data_demo.rds**) for the workflows.
* Run the **STA495MT_JS/STA495MasterThesis/sim_workflow_tikz/flow_para.Rnw** and **STA495MT_JS/STA495MasterThesis/sim_workflow_tikz/flow_design.Rnw** files to generate the workflows.
* Compile (*Build All*) the **STA495MT_JS/STA495MasterThesis/report/report.Rproj** project again to finalize the report.

### What is required to work with these repositories? ###
This thesis was conducted with the R version 4.2.2 (2022-11-10) under Ubuntu Linux 20.04 LTS and **knitr** is used for compiling. Use the following code to install all the required R packages 

* install.packages(c("knitr","tram", "RColorBrewer", "tidyverse", "ggtext","ggside", "colorspace", "mlbench", "xtable","gridExtra", "metR", "tableone", "scatterplot3d", "scales","mvtnorm", "survival"))

The **biostatUZH** and **Collinearity** packages are also needed but are not available from CRAN. They can be installed as:

* install.packages("biostatUZH", repos="http://R-Forge.R-project.org")
* install.packages("devtools") 
* library("devtools")
* devtools::install_github(repo = "jsepin/Collinearity") 

### Who do I talk to if I have questions? ###

* jerome.sepin@uzh.ch
* jerome.sepin@gmail.com
