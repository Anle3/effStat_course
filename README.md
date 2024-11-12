                                             
<img src="readme.PNG" width="200" align="right" />

# Machine Learning in clinical development README

## Introduction

The two demos that we have build aim at providing basic introduction to ML with R for the survival outcome using random forest (demo_ML_survivalRF.html) and the application of virtual twins for exploring treatment effect heterogeneity (demo_VirtualTwins.html). We used simulated data. 

# Prerequisites

All code was generated using Rstudio with R version 4.2.1 (for demo_ML_survivalRF.Rmd) or higher (for demo_VirtualTwins.html).  Rendering the R Markdown document (.Rmd) file using this or newer versions in Rstudio should not be a problem.

# Open the project 

From inside RStudio go to: "File-> Open Project" and select: path_to_ml_course/effStat_course/effStat_course.Rproj

This will open the project of the course

## Use of renv package (specifically for demo_ML_survivalRF)

In order to ensure seamless rendering of the demo_ML_survivalRF.Rmd file you need to ensure you use the same environment as the one we used when we created the file. Thus we have used **renv**, a dependency management package. 
You might need to install packages recorded in the renv lockfile.
You should use the following commands

```
install.packages("renv") #in case you dont have renv installed
renv::restore()
```
[renv tutorial](https://rstudio.github.io/renv/articles/renv.html)

## Pre-processed data

**./Data** folder

data_prep.qs is the actual file containing the pre-processed simulated data that has been generated for both demos. 

## Model location

**./ML** folder

htune_demo.qs is the actual file containing the RF ML model that has been generated for the first demo (demo_ML_survivalRF). 

## Location of the .RMD file

The script needed to run the first demo is:


-demo_ML_survivalRF.Rmd, the actual .Rmd file

-demo_ML_survivalRF.html the output fom rendering the .Rmd

The script needed to run the second demo is:

-demo_VirtualTwins.Rmd, the actual .Rmd file

-demo_VirtualTwins.html the output fom rendering the .Rmd

## Rendering Output

There are two ways to render an R Markdown document. If you are using RStudio, then the “Knit” button (Ctrl+Shift+K) will render the document. Alternatively you can use the run_all.R script.

# Troubleshooting

If any questions feel free to contact:


- Sandra Gonzalez Maldonado: <sandra.gonzalez_maldonado@boehringer-ingelheim.com>

- Eliana Garcia Cossio: <eliana.garciacossio@bayer.com>

- Antigoni Elefsinioti: <antigoni.elefsinioti@bayer.com>


  
