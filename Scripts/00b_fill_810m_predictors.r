
# fill in shorelines on 810m climate data, for use in final maxent models
# matthew kling, january 2016

cdir <- climate_data_dir
outdir <- filled_climate_data_dir

fill_iterations <- 3

f <- list.files(cdir, pattern="810m", full.names=T)

# function fills NA values along shorlines with the average of adjacent (queen) values.
# iterations = number of times to run, filling in one layer of pixels each time.
fill <- function(x, iterations=1){
        for(i in 1:iterations) x <- focal(x, w=matrix(1,nrow=3,ncol=3), function(x)mean(x[!is.na(x)]), NAonly=T)
        return(x)
}

for(file in f){
        r <- readRDS(file)
        r <- fill(r, iterations=fill_iterations)
        saveRDS(r, paste0(outdir, "/", sub(".rdata", paste0("_filled", fill_iterations, "x.rds"), basename(file))))
}



