                                             
<img src="readme.PNG" width="200" align="right" />

# SDI ML course README

## Introduction

This Demo aims at providing basic introduction to ML with R. We used simulated data that we generated internally. 

# Prerequisites

All code was generated using Rstudio with R version 4.0.1.  Rendering the R Markdown document (.Rmd) file using this or newer versions in Rstudio should not be a problem

# Open the project 

From inside RStudio go to: "File-> Open Project" and select: path_to_ml_course/sdi_ml_course/sdi_ml_course.Rproj

This will open the project of the course

## Use of renv package

In order to ensure seamless rendering of the demo .Rmd file you need to ensure you use the same environment as the one we used when we created the file. Thus we have used **renv**, a dependency management package. 
You might need to install packages recorded in the renv lockfile.
You should use the following commands

```
install.packages("renv") #in case you dont have renv installed
renv::restore()
```
[renv tutorial](https://rstudio.github.io/renv/articles/renv.html)

## Data location

**./Data** folder

Demo_data_tte.qs is the actual file containing the simulated data. The other files were generated through knitting RMD and are used to save knitting time

## Location of the .RMD file

The script needed to run the package are in the following folder

**./Demo**

-Demo_simulated_tte.Rmd, the actual .Rmd file

-Demo_simulated_tte.html the output fom rendering the .Rmd

## Rendering Output
There are two ways to render an R Markdown document. If you are using RStudio, then the “Knit” button (Ctrl+Shift+K) will render the document.

# Troubleshooting

If any questions feel free to contact

- Antigoni Elefsinioti: <antigoni.elefsinioti@bayer.com>

- Eliana Garcia Cossio: <eliana.garciacossio@bayer.com>
  