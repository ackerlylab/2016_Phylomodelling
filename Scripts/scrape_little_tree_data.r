
# Download all the Little tree map shapefiles from USGS
# Matthew Kling
# March 2016

library(stringr)

# paths to source website and target folder
url <- "http://esp.cr.usgs.gov/data/little/"
destination <- "C:/Lab_projects/2016_Phylomodelling/Data/Species/Little_trees/"

# get unique urls for every species
html <- paste(readLines(url), collapse="\n")
files <- gregexpr("\\.zip", html)[[1]]
files <- sapply(files, function(x) substr(html, x-8, x+3))

# add latin names extracted from web page
latin <- gregexpr("<I>", html)[[1]] # this works because these are the only italicized elements on the page
latin <- sapply(latin, function(x) substr(html, x+3, x+50))
latin <- as.vector(sapply(latin, function(x) substr(x, 1, regexpr("</I>", x)[[1]]-1)))
files <- cbind(files, latin)
if(all.equal(toupper(substr(files[,1],1,2)), substr(files[,2],1,2)) == F) stop("name matching failed")

# download files
getfile <- function(x){
        x <- sapply(x, str_trim)
        zipdest <- paste0(destination, x[1])
        download.file(paste0(url, x[1]), zipdest)
        unzip(zipdest, exdir=paste0(destination, x[2])) # unzip into latin named folder
        file.remove(zipdest)
}
apply(files, 1, getfile)