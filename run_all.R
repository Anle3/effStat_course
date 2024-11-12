path_pgms<- "~/effStat_course2/"

# first demo
rmarkdown::render(                
  input       = fs::path(path_pgms, "demo_ML_survivalRF.Rmd"),
  output_dir  = path_pgms,
  envir       = new.env()
)

# second demo
rmarkdown::render(                
  input       = fs::path(path_pgms, "demo_VirtualTwins.Rmd"),
  output_dir  = path_pgms,
  envir       = new.env()
)