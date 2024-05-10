
library(phyloseq)
library(ggplot2)
library(ggbeeswarm)
library(rcompanion)
library(FSA)
library(gtools)
# load data
load("results/4.Rarefy_and_subset/all_raref.RData")
load("results/4.Rarefy_and_subset/sal_raref.RData")
load("results/4.Rarefy_and_subset/env_raref.RData")
metadata <- data.frame(sample_data(all_r))
# calculate alpha diversity
env.richness <- estimate_richness(env_r, measures = "Observed")
sal.richness <- estimate_richness(sal_r, measures = "Observed")
# add extra metadata cols
mdcols <- c("Species", "Locality", "Habitat")
for(C in mdcols){
  env.richness[,C] <- metadata[rownames(env.richness),C]
  sal.richness[,C] <- metadata[rownames(sal.richness),C]
}
env.richness$Habitat <- paste(" ", env.richness$Habitat, sep = "")
env.richness$Species <- env.richness$Habitat
## Scheirer-Ray-Hare tests
# env samples: test for Locality and Habitat
# subset data
data.env <- data.frame(Habitat = factor(env.richness$Habitat), 
                       Locality = factor(env.richness$Locality), 
                       Observed = env.richness$Observed)
# SRH test
a <- "Observed"

SRH.env <- scheirerRayHare(formula(paste(a, "~ Locality * Habitat")), data = data.env)
# sal samples: only one species (N. vir) occurs in more than one habitat, so split data by habitat and test each subset for Species & Locality. Also subset data for N. vir and test that subset for Habitat and Locality.
# habitat subsets
data.sal <- list()
habs <- sort(unique(sal.richness$Habitat))
for(h in habs) {
  data.sal[[h]] <- data.frame(Species = factor(sal.richness[which(sal.richness$Habitat == h),"Species"]), 
                              Habitat = factor(sal.richness[which(sal.richness$Habitat == h),"Habitat"]), 
                              Locality = factor(sal.richness[which(sal.richness$Habitat == h),"Locality"]), 
                              Observed = sal.richness[which(sal.richness$Habitat == h),"Observed"])
}
# N. vir subsets
data.sal[["Nvir"]] <- data.frame(Species = factor(sal.richness[which(sal.richness$Species == "N. viridescens"),"Species"]),
                                 Habitat = factor(sal.richness[which(sal.richness$Species == "N. viridescens"),"Habitat"]),
                                 Locality = factor(sal.richness[which(sal.richness$Species == "N. viridescens"),"Locality"]),
                                 Observed = sal.richness[which(sal.richness$Species == "N. viridescens"),"Observed"])
# SRH tests
SRH.sal <- list()
for(h in habs){
  SRH.sal[[h]] <- scheirerRayHare(formula(paste(a, "~ Species * Locality")), data = data.sal[[h]])
}
SRH.sal[["Nvir"]] <- scheirerRayHare(formula(paste(a, "~ Habitat * Locality")), data = data.sal[["Nvir"]])

## Dunn's posthoc tests for significant variables in SRH tests
# environmental samples
DunnTests <- list()
my.vars <- c("Locality", "Habitat")
for(v in my.vars){
  if(SRH.env[v, "p.value"] <= 0.05 & length(unique(data.env[,v])) > 2){
      my.test <- dunnTest(formula(paste(a, "~", v)), data = data.env, method = "bonferroni")$res
      my.test$dataset <- "env"
      my.test$factor <- v 
      DunnTests[[paste0("env.",v)]] <- my.test
  }
}

# sal per habitat
my.vars <- c("Species", "Locality")
for(h in habs){
    for(v in my.vars){
      if(SRH.sal[[h]][v, "p.value"] <= 0.05  & length(unique(data.sal[[h]][,v])) > 2){
        my.test <- dunnTest(formula(paste(a, "~", v)), data = data.sal[[h]], method = "bonferroni")$res
        my.test$dataset <- paste("sal", h, sep = ".")
        my.test$factor <- v 
        DunnTests[[paste("sal", h, v, sep = ".")]] <- my.test
      }
    }
}
# sal N. vir
my.vars <- c("Habitat", "Locality")
  for(v in my.vars){
    if(SRH.sal[["Nvir"]][v, "p.value"] <= 0.05 & length(unique(data.sal[["Nvir"]][,v])) > 2){
    my.test <- dunnTest(formula(paste(a, "~", v)), data = data.sal[["Nvir"]], method = "bonferroni")$res
    my.test$dataset <- "sal.Nvir"
    my.test$factor <- v 
    DunnTests[[paste("sal.Nvir", v, sep = ".")]] <- my.test
    }
}
# combine and recalculate adjusted pvalues
DunnTests.comb <- do.call(rbind, DunnTests)
DunnTests.comb$P.adj <- p.adjust(DunnTests.comb$P.unadj, method = "fdr")
DunnTests.comb$sig <- stars.pval(DunnTests.comb$P.adj)
# save results
dir.create("results/6.Alpha_diversity")
# Dunn tests
write.csv(DunnTests.comb, file = "results/6.Alpha_diversity/DunnTests.csv", row.names = F)
# SRH tests
sink("results/6.Alpha_diversity/SRHTests.txt")
print("Environmental samples")
SRH.env
print("Salamander skin samples")
SRH.sal
sink()
## plot
# combine richnesses
all.richness <- rbind(sal.richness, env.richness)
colnames(all.richness)[1] <- "ASV Richness"
# colours
sppcol  <- c("black","olivedrab","forestgreen","#007054","lightblue","gray","gold","#F2A181","#FFB0C9","darkmagenta")  
names(sppcol) <- sort(unique(metadata$Species))
# function to make italic legend for species
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
# plot ASV richness
p.all.obs <- ggplot(all.richness, aes(x = Habitat, y = `ASV Richness`, color = Species, shape = Locality)) +
  geom_beeswarm(cex = 3, corral = "wrap") + scale_color_manual(values = sppcol, limits = names(sppcol), labels = make.italic(names(sppcol))) +
  theme_test() +
  theme(legend.text.align = 0) + 
  stat_summary(all.richness, fun = mean, geom = "crossbar",mapping = aes(x = Habitat, y = `ASV Richness`, group =1))
p.all.obs

