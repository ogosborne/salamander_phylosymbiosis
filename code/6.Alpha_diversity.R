
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
set.seed(8916075)
env.richness <- estimate_richness(env_r, measures = c("Observed", "Chao1"))
sal.richness <- estimate_richness(sal_r, measures = c("Observed", "Chao1"))
# add extra metadata cols
mdcols <- c("Species", "Locality", "Habitat")
for(C in mdcols){
  env.richness[,C] <- metadata[rownames(env.richness),C]
  sal.richness[,C] <- metadata[rownames(sal.richness),C]
}
env.richness$Habitat <- paste(" ", env.richness$Habitat, sep = "")
env.richness$Species <- env.richness$Habitat
#sal.richness$Habitat <- paste("", sal.richness$Habitat, sep = "")
## Scheirer-Ray-Hare tests
# env samples: test for Locality and Habitat
# subset data
data.env <- data.frame(Habitat = factor(env.richness$Habitat), 
                       Locality = factor(env.richness$Locality), 
                       Observed = env.richness$Observed, 
                       Chao1 = env.richness$Chao1
                       )
# SRH test
a.stats <- c("Observed", "Chao1")
SRH.env <- list()
for(a in a.stats){
  SRH.env[[a]] <- scheirerRayHare(formula(paste(a, "~ Locality * Habitat")), data = data.env)
}
# sal samples: only one species (N. vir) occurs in more than one habitat, so split data by habitat and test each subset for Species & Locality. Also subset data for N. vir and test that subset for Habitat and Locality.
# habitat subsets
data.sal <- list()
habs <- sort(unique(sal.richness$Habitat))
for(h in habs) {
  data.sal[[h]] <- data.frame(Species = factor(sal.richness[which(sal.richness$Habitat == h),"Species"]), 
                              Habitat = factor(sal.richness[which(sal.richness$Habitat == h),"Habitat"]), 
                              Locality = factor(sal.richness[which(sal.richness$Habitat == h),"Locality"]), 
                              Observed = sal.richness[which(sal.richness$Habitat == h),"Observed"], 
                              Chao1 = sal.richness[which(sal.richness$Habitat == h),"Chao1"])
}
# N. vir subsets
data.sal[["Nvir"]] <- data.frame(Species = factor(sal.richness[which(sal.richness$Species == "N. viridescens"),"Species"]),
                                 Habitat = factor(sal.richness[which(sal.richness$Species == "N. viridescens"),"Habitat"]),
                                 Locality = factor(sal.richness[which(sal.richness$Species == "N. viridescens"),"Locality"]),
                                 Observed = sal.richness[which(sal.richness$Species == "N. viridescens"),"Observed"], 
                                 Chao1 = sal.richness[which(sal.richness$Species == "N. viridescens"),"Chao1"])
# SRH tests
SRH.sal <- list()
for(h in habs){
  SRH.sal[[h]] <- list()
  for(a in a.stats){
    SRH.sal[[h]][[a]] <- scheirerRayHare(formula(paste(a, "~ Species * Locality")), data = data.sal[[h]])
  }
}
SRH.sal[["Nvir"]] <-  list()
for(a in a.stats){
  SRH.sal[["Nvir"]][[a]] <- scheirerRayHare(formula(paste(a, "~ Habitat * Locality")), data = data.sal[["Nvir"]])
}
## Dunn's posthoc tests for significant variables in SRH tests
DunnTests <- list()
for(a in a.stats){
  DunnTests[[a]] <- list()
}

# environmental samples
my.vars <- c("Locality", "Habitat")
for(a in a.stats){
  for(v in my.vars){
    if(SRH.env[[a]][v, "p.value"] <= 0.05 & length(unique(data.env[,v])) > 2){
      my.test <- dunnTest(formula(paste(a, "~", v)), data = data.env, method = "bonferroni")$res
      my.test$dataset <- "env"
      my.test$factor <- v 
      DunnTests[[a]][[paste0("env.",v)]] <- my.test
    }
  }
}

# sal per habitat
my.vars <- c("Species", "Locality")
for(a in a.stats){
  for(h in habs){
    for(v in my.vars){
      if(SRH.sal[[h]][[a]][v, "p.value"] <= 0.05  & length(unique(data.sal[[h]][,v])) > 2){
        my.test <- dunnTest(formula(paste(a, "~", v)), data = data.sal[[h]], method = "bonferroni")$res
        my.test$dataset <- paste("sal", h, sep = ".")
        my.test$factor <- v 
        DunnTests[[a]][[paste("sal", h, v, sep = ".")]] <- my.test
      }
    }
  }
}
# sal N. vir
my.vars <- c("Habitat", "Locality")
for(a in a.stats){
  for(v in my.vars){
    if(SRH.sal[["Nvir"]][[a]][v, "p.value"] <= 0.05 & length(unique(data.sal[["Nvir"]][,v])) > 2){
    my.test <- dunnTest(formula(paste(a, "~", v)), data = data.sal[["Nvir"]], method = "bonferroni")$res
    my.test$dataset <- "sal.Nvir"
    my.test$factor <- v 
    DunnTests[[a]][[paste("sal.Nvir", v, sep = ".")]] <- my.test
    }
  }
}
# combine and recalculate adjusted pvalues
for(a in a.stats){
  DunnTests[[a]] <- do.call(rbind, DunnTests[[a]])
  DunnTests[[a]]$P.adj <- p.adjust(DunnTests[[a]]$P.unadj, method = "fdr")
  DunnTests[[a]]$sig <- stars.pval(DunnTests[[a]]$P.adj)
}
DunnTests.comb <- merge(DunnTests$Observed, DunnTests$Chao1, by = c("dataset", "factor", "Comparison"), all = T, suffixes = c(".Observed", ".Chao1"))
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
  theme(legend.text.align = 0)
p.all.obs
# plot Chao1
p.all.chao <- ggplot(all.richness, aes(x = Habitat, y = Chao1, color = Species, shape = Locality)) +
  geom_beeswarm(cex = 3, corral = "wrap") + scale_color_manual(values = sppcol, limits = names(sppcol), labels = make.italic(names(sppcol))) + 
  theme_test() +
  theme(legend.text.align = 0)
p.all.chao

