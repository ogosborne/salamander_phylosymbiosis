library(phyloseq)
library(ggplot2)
library(dplyr)
library(ACAT)
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
#env_abd_adiv <- estimate_richness(env_r_antibd, measures = c("Observed"))
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
  dplyr::group_by(Species) %>%
  dplyr::summarize(Bd_load = mean(Bd_load, na.rm = T), 
            all_adiv = mean(all_adiv), 
            abd_adiv = mean(abd_adiv), 
            abd_relabund = mean(abd_relabund)) %>%
  as.data.frame()
# add prevalence
prev <- as.matrix(table(df_sal[,c("Bd", "Species")]))
prev <- sweep(prev,2,colSums(prev),`/`)["yes",]
#add
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

### measure phylosignal accounting for intra-specific variation using bootstrapping (random sampling with replacement)

### Functions for bootstraping
# get a list containing a vector of sample names for each species
get_samplist <- function(metadata, samp.col, grp.col){
  grps <- sort(unique(metadata[, grp.col]))
  out <- list()
  for(g in grps){
    out[[g]] <- metadata[ which(metadata[, grp.col] == g), samp.col]
  }
  out
}
# get a bootstrap sample with the correct number of individuals per species
bootstrap_per_grp <- function(samp.list){
  # samp.list is a named vector containing one random representative sample name for each group, named by group.
  out <- list()
  for(g in names(samp.list)){
    n <- length(samp.list[[g]])
    out[[g]] <- sample(samp.list[[g]], n, replace = T)
  }
  out <- do.call(c, out)
  out <- unname(out)
  out
}

# run 1000 replicates of multiPhylosignal with bootstrap replicates
df_sal$Samp <- rownames(df_sal)
samplist <- get_samplist(df_sal, "Samp", "Species")
bs.mps <- list()
set.seed(51089)
for(n in 1:1000){
  # sample one ind per species
  my.bs.df <- df_sal[bootstrap_per_grp(samplist),c("Species","Bd","Bd_load","all_adiv","abd_adiv","abd_relabund")]
  # recalculate mean parameters
  my.params <- my.bs.df %>%
    dplyr::group_by(Species) %>%
    dplyr::summarize(Bd_load = mean(Bd_load, na.rm = T), 
              all_adiv = mean(all_adiv), 
              abd_adiv = mean(abd_adiv), 
              abd_relabund = mean(abd_relabund)) %>%
    as.data.frame()
  # add prevalence
  my.prev <- as.matrix(table(my.bs.df[,c("Bd", "Species")]))
  my.prev <- sweep(my.prev,2,colSums(my.prev),`/`)["yes",]
  my.params$Bd_prev <- my.prev
  # reorder paramaters
  row.names(my.params) <- my.params$Species
  my.params$Species = NULL
  order <- match(host.phylo$tip.label, rownames(my.params))
  my.params <- my.params[order,]
  # Run multiPhylosignal
  bs.mps[[n]] <- multiPhylosignal(my.params, host.phylo, reps = 10000)
}
# get all p values 
allP.bs.mps <- data.frame(matrix(0,ncol = 5, nrow = 1000, dimnames = list(NULL, c("Bd_load","Bd_prev","all_adiv","abd_adiv","abd_relabund")))) 
for(n in 1:1000){
  myres <- bs.mps[[n]]
  for(t in c("Bd_load","Bd_prev","all_adiv","abd_adiv","abd_relabund")){
    allP.bs.mps[n,t] <- myres[t,"PIC.variance.P"] 
  }
}
# pool P values with Cauchy method
mps$BS.cauchy.P <- apply(allP.bs.mps, MARGIN=2, FUN=ACAT)[rownames(mps)]
# save
dir.create("results/12.Bd_phylosignal")
write.csv(mps, file = "results/12.Bd_phylosignal/multiphylosignal.csv")
save(bs.mps, file = "results/12.Bd_phylosignal/bootstrap.multiphylosignal.RData")
