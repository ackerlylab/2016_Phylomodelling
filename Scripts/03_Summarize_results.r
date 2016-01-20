# Script to summarize results across all species and save csv files with maxent threshold index for each occurrence record

# Clear workspace
rm(list=ls())

#-----------------#
# Set directories #
#-----------------#
Computer <- 'Ackerly'

if(Computer == 'Ackerly') {
  wdir <- 'E:/Phylo_modelling/'
  cdir <- 'E:/BCM/CA_2014/Summary/HST/Normals_30years/'
} else if (Computer == 'HP') {
  wdir <- 'D:/Phylo_modelling/'
  cdir <- 'D:/BCM/CA_2014/Summary/HST/Normals_30years/'
}

spdir <- paste(wdir, 'Data/Species/Processed/', sep='')
odir <- paste(wdir, 'Output/', sep='')

# What species?
allSpecies <- sort(readRDS(paste(spdir, '0_Species_list.rdata', sep='')))
spp_dirs <- list.dirs("E:/Phylo_modelling/Output/Maxent/V2", recursive=F)

res <- list()
res2 <- list()
#for(i in 1:10) {
for(i in 1:length(allSpecies)) {
  x <- allSpecies[i]
  if (!file.exists(paste0("E:/Phylo_modelling/Output/Maxent/V2/", x, "/Suitability_data.rdata"))) {
    next()
  }
  d <- readRDS(paste0("E:/Phylo_modelling/Output/Maxent/V2/", x, "/Suitability_data.rdata"))
  f <- d$occur
  f$prop_max.suit.occ[f$prop_max.suit.occ<10e-4] <- 10e-4 # put a lower limit on the ratios
  f$outlier <- -log10(f$prop_max.suit.occ)
  res[[i]] <- f[,c('id','current_name_binomial','outlier')]
  res2[[i]] <- d$summary
}
result <- do.call("rbind", res)

result1 <- result[rev(order(result$outlier)),]
result2 <- do.call("rbind",res2)

write.csv(result1, paste0(odir,'Maxent_outlier_indices.csv'))
write.csv(result2, paste0(odir,'Maxent_summary_stats.csv'))


library(ggplot2)
p <- ggplot(result2, aes(records.used, auc.train)) + 
  geom_point(color="lightblue") +
  stat_density2d(color="black", size=1) +
  scale_x_log10(breaks=c(1, 10, 100, 1000)) +
  theme_minimal() +
  labs(x="number of records", y="training AUC", 
       title=paste0("Maxent training AUC vs. number of records\n(n = ", nrow(result2), " plant species)"))
ggsave("E:/Phylo_modelling/Output/Charts/auc_vs_records.png", p, width=8, height=6)

