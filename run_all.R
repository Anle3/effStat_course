path_pgms<- "~/effStat_course2/"
rmarkdown::render(                
  input       = fs::path(path_pgms, "demo.Rmd"),
  output_dir  = path_pgms,
  envir       = new.env()
)