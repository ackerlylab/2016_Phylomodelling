
#### attempt to quantify, understand, visualize niche model's marginal improvement over geodistance null


# set directories
cdir <- filled_climate_data_dir
spdir <- occurrence_data_dir_processed
odir <- "C:/Lab_projects/2016_Phylomodelling/Output/Niche"
species_list <- paste0(spdir, 'combined/0_Species_list_v2.rdata')
maxent_background <- paste(spdir, 'combined/10000_CA_plants_bg_810m_occbias.rdata', sep='')

# load libraries
options(java.parameters = "-Xmx1g" )
require(dismo)
require(rJava)
require(raster)
library(rgeos)
library(rgdal)
extract <- raster::extract

# load climate data
files <- list.files(path=cdir, pattern='810m_filled3x', full.names=FALSE)
climnames <- sort(c("cwd","djf","jja","ppt")) #use sort to ensure correct names used for predictors
predictors <- stack(lapply(files,function(x) readRDS(paste(cdir,x,sep="/"))))
names(predictors) <-  climnames 

# load occurrences
allocc <- list.files("C:/Lab_projects/2016_Phylomodelling/Data/Species/processed3/atomic", full.names=T)
allocc <- lapply(allocc, readRDS)
allocc <- do.call("rbind", allocc)
allocc <- allocc[!is.na(allocc$longitude + allocc$latitude),]
coordinates(allocc) <- c("longitude", "latitude")
projection(allocc) <- '+proj=longlat +ellps=WGS84'
allocc <- spTransform(allocc, crs(predictors))

# exclude occurrences with no climate data
allocc <- allocc[!is.na(raster::extract(predictors[[1]], allocc)),]
species <- unique(allocc$current_name_binomial)

# raster of unique pixel IDs
pixels <- predictors[[1]]
pixels[] <- 1:ncell(pixels)



process <- function(spp){
        
        require(raster)
        require(dismo)
        require(rJava)
        require(rgeos)
        require(rgdal)
        
        # maxent background
        bg_all <- readRDS(maxent_background)
        
        # unique occurrences of species
        occ <- allocc[allocc$current_name_binomial==spp,]
        occ <- occ[sample(nrow(occ), nrow(occ), replace=F),] # scramble before getting unique, to prevent systematic within-pixel directional bias when creating xval folds
        occ <- occ[!duplicated(extract(pixels, occ)),] # get unique
        
        if(nrow(occ)<20) return("insufficient points")
        
        # select background points unique from each other and from occurrences
        bg <- bg_all
        bg <- bg[!duplicated(extract(pixels, bg)),]
        bg <- bg[!extract(pixels, bg) %in% extract(pixels, occ),]
        
        # assign presences to xval folds
        occ$fold <- sample(4, nrow(occ), replace=T)
        bg$fold <- sample(4, nrow(bg), replace=T)
        
        s <- data.frame()
        for(fold in 1:4){
                
                # training points
                train_pres <- occ[occ$fold != fold,]
                train_bg <- bg[bg$fold != fold,]
                
                # evaluation points
                test_pres <- occ[occ$fold == fold,]
                test_bg <- bg[bg$fold == fold,]
                
                # subset eval points via pairwise distance sampling
                selected <- pwdSample(fixed=test_pres, sample=test_bg, reference=train_pres, tr=.2, lonlat=F)
                test_bg <- test_bg[as.vector(na.omit(selected)),]
                test_pres <- test_pres[!is.na(as.vector(selected)),]
                
                # fit maxent model and make prediction
                mxArgs <- c("-a", "-z", "outputformat=raw", "maximumbackground=10000", "nothreshold", "nohinge")
                mx <- maxent(predictors, p=train_pres, args=mxArgs, a=train_bg)
                pred <- predict(mx, predictors)
                
                # evaluate
                eval <- evaluate(test_pres, test_bg, mx, predictors)
                
                # extract data for test points
                rstrs <- stack(pred, predictors)
                names(rstrs)[1] <- "maxent"
                pres <- as.data.frame(extract(rstrs, test_pres))
                abs <- as.data.frame(extract(rstrs, test_bg))
                pres$presence="presence"
                abs$presence="absence"
                pres <- cbind(pres, coordinates(test_pres)) #############
                abs <- cbind(abs, coordinates(test_bg))
                
                # compile
                fd <- rbind(pres, abs)
                fd$fold <- fold
                fd$species <- spp
                fd$auc <- eval@auc
                fd$nocc <- nrow(occ)
                s <- rbind(s, fd)
        }
        return(s)
}


library(doParallel)
cl <- makeCluster(7)
registerDoParallel(cl)
results <- foreach(spp=species,
                   .packages=c("raster", "dismo", "rJava", "rgeos", "rgdal")) %dopar% {
                           return(try(process(spp)))
                   }
stopCluster(cl)

# compile and save results
s <- do.call("rbind", results[sapply(results, class)=="data.frame"])
saveRDS(s, paste0(odir, "/pwds_climate_results.rds"))


stop("data prep complete")


##### viz & analysis #####

# currently this is only CA endemics that occur in 20+ cells

library(tidyr)
library(dplyr)
library(broom)
library(ggplot2)
library(gclus)


# restrict to cali endemics
e <- read.csv("C:/Lab_projects/2016_Phylomodelling/Data/Endemism/Californian_endemic_binomials.txt", 
              stringsAsFactors=F, header=F)[,1]
s <- readRDS(paste0(odir, "/pwds_climate_results.rds"))
s <- s[s$species %in% e,]
s$presence <- factor(s$presence)

# distribution of AUC scores across taxa
u <- s %>%
        group_by(species) %>%
        summarize(auc=mean(auc))
p <- ggplot(u, aes(auc)) + 
        geom_density(fill="darkred", alpha=.5, color=NA) +
        geom_vline(xintercept=c(.5,1)) +
        labs(title="PWDS AUC values across CA endemic species") +
        theme_minimal()
ggsave(paste0(odir, "/auc_density.png"), p, width=8, height=6)


# plot pres-abs in multivariate climate space
# compare this to the base maxent model results with no pwds
# 4 90% ellipses: maxent pred vs bg; pwds pres vs abs

aucs <- s %>%
        group_by(species, fold) %>%
        summarize(auc=mean(auc)) %>%
        group_by(species) %>%
        summarize(auc=mean(auc))

myspp <- sample(unique(s$species), 25)
p <- ggplot() + 
        stat_ellipse(data=select(s, -species), aes(djf, log10(ppt)), 
                     level=.999, geom="polygon", fill="gray", color=NA) +
        stat_ellipse(data=s[s$species %in% myspp,], aes(djf, log10(ppt), fill=presence), 
                     level=.95, size=1, alpha=.5, color=NA, geom="polygon") +
        geom_text(data=aucs[aucs$species %in% myspp,], aes(x=mean(s$djf),y=mean(log10(s$ppt)),label=round(auc, 2),size=auc)) +
        theme_minimal() +
        facet_wrap(~species) +
        labs(title="PWDS: 95% ellipses for 25 random CA endemics")
ggsave(paste0(odir, "/ellipses_ppt_djf_25.png"), p, width=12, height=9)



##################
# VARIABLE IMPORTANCE
##################

# WILCOXON
v <- s %>%
        mutate(ppt = log10(ppt)) %>%
        gather(variable, value, cwd:ppt) %>%
        group_by(species, variable) %>%
        summarize(wp = wilcox.test(value[presence=="presence"], value[presence=="absence"])$p.value)
p <- ggplot(v, aes(variable, wp, fill=variable)) + 
        geom_boxplot(alpha=.5) +
        geom_hline(yintercept=.05, color="red") +
        theme_minimal() +
        scale_y_sqrt(breaks=c(.01,.1,.5,1)) +
        labs(y="wilcoxon p-value",
             title="WILCOXON P-VALUES") +
        theme(legend.position="none")
ggsave(paste0(odir, "/variables_boxplot_wilcoxon_p.png"), p, width=8, height=6)

# UNIVARIATE REGRESSION
v <- s %>%
        mutate(ppt = log10(ppt)) %>%
        mutate_each(funs(scale), cwd:ppt) %>%
        gather(variable, value, cwd:ppt) %>%
        group_by(species, variable) %>%
        do(tidy(glm(presence ~ value, data=., family="binomial")))

p <- ggplot(v[v$term != "(Intercept)",], aes(variable, p.value, fill=variable)) + 
        geom_boxplot(alpha=.5) +
        geom_hline(yintercept=.05, color="red") +
        theme_minimal() +
        scale_y_sqrt(breaks=c(.01,.1,.5,1)) +
        labs(y="univariate regression p-value",
             title="UNIVARIATE REGRESSION P-VALUES") +
        theme(legend.position="none")
ggsave(paste0(odir, "/variables_boxplot_unireg_p.png"), p, width=8, height=6)

p <- ggplot(v[v$term != "(Intercept)",], aes(variable, abs(estimate), fill=variable)) + 
        geom_boxplot(alpha=.5) +
        theme_minimal() +
        scale_y_log10() +
        labs(y="univariate regression p-value",
             title="UNIVARIATE LOGISTIC REGRESSION COEFFIEINTS") +
        theme(legend.position="none")
ggsave(paste0(odir, "/variables_boxplot_unireg_slope.png"), p, width=8, height=6)


########## MULTIPLE REGRESSION ###########
v <- s %>%
        mutate(ppt = log10(ppt)) %>%
        mutate_each(funs(scale), cwd:ppt) %>%
        group_by(species) %>%
        do(tidy(glm(presence ~ cwd + djf + jja + ppt, data=., family="binomial")))

p <- ggplot(v[v$term != "(Intercept)",], aes(term, p.value, fill=term)) + 
        geom_boxplot(alpha=.5) +
        geom_hline(yintercept=.05, color="red") +
        theme_minimal() +
        scale_y_sqrt(breaks=c(.01,.1,.5,1)) +
        labs(y="multiple regression p-value",
             title="MULTIPLE LOGISTIC REGRESSION P-VALUES") +
        theme(legend.position="none")
ggsave(paste0(odir, "/variables_boxplot_multreg_p.png"), p, width=8, height=6)

p <- ggplot(v[v$term != "(Intercept)",], aes(term, abs(estimate), fill=term)) + 
        geom_boxplot(alpha=.5) +
        theme_minimal() +
        scale_y_log10() +
        labs(y="absolute value of multiple regression coefficient",
             title="MULTIPLE LOGISTIC REGRESSION COEFFICIENTS") +
        theme(legend.position="none")
ggsave(paste0(odir, "/variables_boxplot_multreg_slope.png"), p, width=8, height=6)


# 2d density plot -- relationships between coefficients and p-values

# map -- spatial variation in variable importance

# scatterplot matrix -- relationships among coefficients
w <- v %>%
        select(species, term, estimate) %>%
        spread(term, estimate) %>%
        select(cwd:ppt) %>%
        as.matrix()
w[abs(w)>5] <- NA
png(paste0(odir, "/variables_scatterplot_multreg_slope.png"), width=1000, height=1000)
pairs(w, pch=16, lwd=2,
      upper.panel=panel.smooth, lower.panel=panel.smooth, span=.2, iter=500,
      main="MULTIPLE LOGISTIC REGRESSION COEFFICIENTS")
dev.off()


##################


d <- data.frame()
for(spp in sample(unique(s$species), 25)){
        sd <- s[s$species==spp,]
        sd$ppt <- log10(sd$ppt)
        pc <- prcomp(sd[,c("cwd", "djf", "jja", "ppt")], center=T, scale=T)
        sd <- cbind(sd, pc$x)
        d <- rbind(d, sd)
}




p <- ggplot() + 
        stat_ellipse(data=select(s, -species), aes(PC1, PC2), level=.99, geom="polygon", fill="gray", color=NA) +
        stat_ellipse(data=d, aes(PC1, PC2, fill=presence), level=.90, size=1, alpha=.5, color=NA, geom="polygon") +
        geom_text(data=aucs, aes(x=0,y=0,label=round(auc, 2),size=auc)) +
        theme_minimal() +
        facet_wrap(~species)
p



p <- ggplot(d) + 
        geom_point(aes(PC1, PC2, alpha=maxent), size=2) +
        stat_density2d(aes(PC1, PC2, color=presence), size=1) +
        theme_minimal() +
        facet_wrap(~species, scales="free")
p

