library(phyloseq)
library(ggplot2)
library(dplyr)
library(ape)
library(picante)
### prep data
# load data
load("results/4.Rarefy_and_subset/sal_raref.Rdata")
load("results/4.Rarefy_and_subset/env_raref.Rdata")
# load antiBD 
AntiBd <- read.csv("results/5.Antifungal_DB/blast.filt.uniq.csv", header = T)
# make antiBd only datasets
sal_r_antibd <-  subset_taxa(sal_r, ASV %in% AntiBd$qseqid)
env_r_antibd <-  subset_taxa(env_r, ASV %in% AntiBd$qseqid)
# Calc alpha diversity overall and for anti-Bd dataset in sal and env
sal_all_adiv <- estimate_richness(sal_r, measures = c("Observed"))
sal_abd_adiv <- estimate_richness(sal_r_antibd, measures = c("Observed"))
env_all_adiv <- estimate_richness(env_r, measures = c("Observed"))
env_abd_adiv <- estimate_richness(env_r_antibd, measures = c("Observed"))
# Relative abundance of antiBd
sal_abd_relabund <- sample_sums(sal_r_antibd) / sample_sums(sal_r)
env_abd_relabund <- sample_sums(env_r_antibd) / sample_sums(env_r)
# Get dataframes
df_sal <- as(sample_data(sal_r), "data.frame")
df_env <- as(sample_data(env_r), "data.frame")
# add div + RA to dataframe
df_sal$all_adiv <- sal_all_adiv$Observed
df_sal$abd_adiv <- sal_abd_adiv$Observed
df_sal$abd_relabund <- sal_abd_relabund
df_env$all_adiv <- env_all_adiv$Observed
df_env$abd_adiv <- env_abd_adiv$Observed
df_env$abd_relabund <- env_abd_relabund
### Correlations between richness and anti-Bd richness for sals and environment
sal_rich_lm <- lm(df_sal$abd_adiv ~ df_sal$all_adiv)
summary(sal_rich_lm)
env_rich_lm <- lm(df_env$abd_adiv ~ df_env$all_adiv)
summary(env_rich_lm)
### Prep data for phylosignal analysis
# get per species
params <- df_sal %>%
  group_by(Species) %>%
  summarize(Bd_load = mean(Bd_load, na.rm = T), 
            all_adiv = mean(all_adiv), 
            abd_adiv = mean(abd_adiv), 
            abd_relabund = mean(abd_relabund)) %>%
  as.data.frame()
# add prevalence
prev <- as.matrix(table(df_sal[,c("Bd", "Species")]))
prev <- sweep(prev,2,colSums(prev),`/`)["yes",]
params$Bd_prev <- prev
# load tree
host.phylo <- read.tree("data/timetree_31JUL23/host.timetree.31JUL23_short.nwk")
host.phylo$tip.label <- gsub("_", " ", host.phylo$tip.label)
# reorder params
row.names(params) <- params$Species
params$Species = NULL
order <- match(host.phylo$tip.label, rownames(params))
params <- params[order,]
### Run multiPhylosignal to measure any phyolgenetic signal in the parameters
set.seed(1212)
mps <- multiPhylosignal(params, host.phylo, reps = 10000)
mps
dir.create("results/12.Bd_phylosignal")
write.csv(mps, file = "results/12.Bd_phylosignal/multiphylosignal.csv")

