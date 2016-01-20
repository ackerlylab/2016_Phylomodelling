


infiles <- list.files("E:/Phylo_modelling/Output/Maxent/V2", pattern="outlier_map_v3.png", recursive=T, full.names=T)

outdir <- "E:/Phylo_modelling/Output/Charts/Mapped_outliers"

for(file in infiles){
  dir <- dirname(file)
  species <- substr(dirname(file), tail(gregexpr("/", dir)[[1]], 1)+1, 1000)
  newfile <- paste0(outdir, "/", species, ".png")
  if(file.exists(newfile)) next()
  #writeLines(paste("moved", species))
  file.copy(file, newfile)
}
