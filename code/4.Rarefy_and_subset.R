library(ape)
library(phyloseq)
library(Biostrings)
library(vegan)
library(tidyverse)
library(gridExtra)

load("results/2.Filter_data/filtered.RData")
# import trees and add to phyloseq objects, multi2di resolves any polytomies
phy_tree(all_f) <- multi2di(read_tree("results/3.ASV_tree/all_f.nwk"))
# rarefaction curve
otu_mat <- as(otu_table(all_f), "matrix")
rc <- rarecurve(t(otu_mat), step=100, tidy = T)
colnames(rc) <- c("Sample", "Depth", "N ASVs")
md <- data.frame(sample_data(all_f))
# split sal and env samples
# sal
sal.samples <- md[which(md$Sample_Type == "Salamander"),]
spp.map <- sal.samples$Species
names(spp.map) <- rownames(sal.samples)
rc.sal <- rc[which(rc$Sample %in% rownames(sal.samples)),]
rc.sal$Species <- spp.map[as.character(rc.sal$Sample)]
# env
env.samples <- md[which(md$Sample_Type == "Environment"),]
hab.map <- env.samples$Habitat
names(hab.map) <- rownames(env.samples)
rc.env <- rc[which(rc$Sample %in% rownames(env.samples)),]
rc.env$Habitat <- hab.map[as.character(rc.env$Sample)]
# N reads per sample
ss.df.sal <- sample_sums(all_f)[rownames(sal.samples)]
ss.df.sal <- data.frame(Sample = names(ss.df.sal), N.reads = ss.df.sal)
ss.df.env <- sample_sums(all_f)[rownames(env.samples)]
ss.df.env <- data.frame(Sample = names(ss.df.env), N.reads = ss.df.env)
# col
sppcol  <- c("black","olivedrab","forestgreen","#007054","lightblue","gray","gold","#F2A181","#FFB0C9","darkmagenta")  
names(sppcol) <- sort(unique(sal.samples$Species))
habcol <- c("#009E73", "#0072B2", "#56B4E9")
names(habcol) <- sort(unique(env.samples$Habitat))
# plot
sp.rp <- ggplot(data = rc.sal, aes(x = Depth, y = `N ASVs`, group = Sample, color = Species)) + 
  geom_line() +  
  geom_vline(xintercept=2945, color = "black", lty = 2) +
  scale_color_manual(values = sppcol)  + 
  theme_classic() +
  ggtitle("B")
sp.hist <- ggplot(data = ss.df.sal, aes(x = N.reads)) + 
  geom_histogram() + 
  ylab("N samples") + 
  xlab("N reads") +
  geom_vline(xintercept=2945, color = "black", lty = 2) +
  theme_classic() +
  ggtitle("A")
en.rp <- ggplot(data = rc.env, aes(x = Depth, y = `N ASVs`, group = Sample, color = Habitat)) + 
  geom_line() +  
  geom_vline(xintercept=2945, color = "black", lty = 2) +
  scale_color_manual(values = habcol)  + 
  theme_classic() +
  ggtitle("D")
en.hist <- ggplot(data = ss.df.env, aes(x = N.reads)) + 
  geom_histogram() + 
  ylab("N samples") + 
  xlab("N reads") +
  geom_vline(xintercept=2945, color = "black", lty = 2) +
  theme_classic() +
  ggtitle("C")
grid.arrange(sp.hist, sp.rp, en.hist, en.rp, 
             layout_matrix = matrix(c(1,2,2,3,4,4), ncol=3, byrow = T))
# remove low count samples
all_f <- prune_samples(sample_sums(all_f ) >= 2945, all_f )
# Rarefy
all_r <- rarefy_even_depth(all_f, replace = FALSE, rngseed = 40191)
# env only
env_r <- subset_samples(all_r, Sample_Type == "Environment")
env_r <- filter_taxa(env_r, function(x) sum(x) != 0, TRUE)
env_f <- subset_samples(all_f, Sample_Type == "Environment")
env_f <- filter_taxa(env_f, function(x) sum(x) != 0, TRUE)
# salamanders only
sal_r <- subset_samples(all_r, Sample_Type == "Salamander")
sal_r <- filter_taxa(sal_r, function(x) sum(x) != 0, TRUE)
sal_f <- subset_samples(all_f, Sample_Type == "Salamander")
sal_f <- filter_taxa(sal_f, function(x) sum(x) != 0, TRUE)
# save
dir.create("results/4.Rarefy_and_subset")
save(all_f, file = "results/4.Rarefy_and_subset/all_unrar.RData")
save(all_r, file = "results/4.Rarefy_and_subset/all_raref.RData")
save(sal_f, file = "results/4.Rarefy_and_subset/sal_unrar.RData")
save(sal_r, file = "results/4.Rarefy_and_subset/sal_raref.RData")
save(env_f, file = "results/4.Rarefy_and_subset/env_unrar.RData")
save(env_r, file = "results/4.Rarefy_and_subset/env_raref.RData")
# save sequences to blast to anti-fungal DB
dir.create("results/5.Antifungal_DB")
writeXStringSet(phyloseq::refseq(all_f), "results/5.Antifungal_DB/all_f.fasta")
