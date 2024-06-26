# libraries
library(phyloseq) 
library(phytools) ## for cospeciation(), cophylo()
library(rlist) ## for list.filter() 
library(vegan) # for mantel
library(ape) # for parafit() & general phylogeny tools
library(ecodist) # for MRM()
library(dplyr) # for %>% 
library(Biostrings) # for writeXStringSet()
library(usedist) # for dist_subset
library(gtools) # for stars.pvals()
source("code/phylosymbiosis_funcs.R") # custom functions

######## 1. DATA PREPARATION ########
## load phyloseq objects for sample data
load("results/4.Rarefy_and_subset/sal_raref.RData") 
load("results/4.Rarefy_and_subset/env_raref.RData") # for sample data only
load("results/4.Rarefy_and_subset/all_unrar.RData") # to get global ASVs
## add combined locality-habitat variable
sample_data(sal_r)$Loc_Hab <- gsub(" ", ".", paste(sample_data(sal_r)$Locality, sample_data(sal_r)$Habitat, sep = "_"))
sample_data(sal_r)$sample <- sample_names(sal_r)
sample_data(env_r)$Loc_Hab <- gsub(" ", ".", paste(sample_data(env_r)$Locality, sample_data(env_r)$Habitat, sep = "_"))
sample_data(all_f)$Loc_Hab <- gsub(" ", ".", paste(sample_data(all_f)$Locality, sample_data(all_f)$Habitat, sep = "_"))
## extract metadata
sal.metadata <- data.frame(sample_data(sal_r))
env.metadata <- data.frame(sample_data(env_r))
## load microbiome distance
load("results/7.Beta_diversity/dismats.RData")
## load host phylogeny
host.phylo <- read.tree("data/timetree_31JUL23/host.timetree.31JUL23_short.nwk")
host.phylo$tip.label <- gsub("_", " ", host.phylo$tip.label)
## load env distances 
load("results/8.Environmental_distances/env.dists.RData")
## make output dir
dir.create("results/9.Phylosymbiosis/")

######## 2. IS THERE A SIGNAL OF PHYLOSYMBIOSIS? ########
## names of the 4 beta-diversity statistics calculated
dists <- c("ja", "bc", "uu", "wu")
## test at individual level with mantel test
set.seed(61511023)
mantel.ind.results <- list()
for(d in dists){
  mantel.ind.results[[d]] <- vegan::mantel(xdis = covar_dists$host_dist, ydis = dismats[[paste0("sal.", d)]], permutations = 10000)
}
mantel.ind.summary <- sum.mantel(mantel.ind.results)
write.csv(mantel.ind.summary, file = "results/9.Phylosymbiosis/ind.mantel.csv", row.names = F)
## produce median salamander skin microbiome distance matrices  (one for each distance measure) for mantel tests
med.sal.dist.per.spp <- list() 
for(d in dists){
  mymat <- summarise_mat(x = as.matrix(dismats[[paste0("sal.", d)]]), metadata = sal.metadata, col = "Species", fun = "median", diagzero = TRUE)
  mymat <- as.dist(mymat)
  med.sal.dist.per.spp[[d]] <- mymat
}
rm(list = c("d", "mymat"))
## produce median salamander skin microbiome trees (one for each distance measure) for cospeciation tests
spp.clust <- list()
for(d in dists){
  tree <- ape::nj(med.sal.dist.per.spp[[d]])
  anc <- mrca(tree)["A. jeffersonianum", "N. viridescens"]
  tree <- root(tree, node = anc)
  spp.clust[[d]] <- tree
}
## test for significant similarity between mean microbiome tree and host phylogeny (phylosymbiosis) with:
# mantel test
set.seed(63737408)
mantel.spp.results <- list()
for(d in dists){
  mantel.spp.results[[d]] <- vegan::mantel(xdis = med.sal.dist.per.spp[[d]], ydis = covar_dists_spp$host_dist, permutations = 10000)
}
mantel.spp.summary <- sum.mantel(mantel.spp.results)
write.csv(mantel.spp.summary, file = "results/9.Phylosymbiosis/spp.mantel.csv", row.names = F)
# cospeciation test
set.seed(97300012)
cospec.results <- list()
for(d in dists){
  cospec.results[[d]] <-  cospeciation(host.phylo, spp.clust[[d]], distance = "RF", method = "permutation", assoc = cbind(host.phylo$tip.label, host.phylo$tip.label), nsim = 10000)
}
cospec.summary <- sum.cospec(cospec.results)
write.csv(cospec.summary, file = "results/9.Phylosymbiosis/RF.perm.test.csv", row.names = F)
## plot cophylo
sppcol <- c("black","olivedrab","forestgreen","#007054","lightblue","gray","gold","#F2A181","#FFB0C9","darkmagenta")
names(sppcol) <- sort(unique(sal.metadata$Species))
for(d in dists){
  plot.col.cophylo(t1 = host.phylo, t2 = spp.clust[[d]], col = sppcol)
}

######## 3. DO ENVIRONMENTAL COVARIATES CONTRIBUTE TO PHYLOSYMBIOSIS? ########
# test the contributions of host phylogenetic distance, environmental microbiome distance, and geographic distance to patterns of phylosymbiosis
# make covars into lower dists
## MRM
# standardise distance matrices
dismats_mrm <- dismats[paste0("sal.", dists)]
covar_dists_d <- lapply(covar_dists, function(x) as.dist((x - mean(x)) / sd(x)))
dismats_mrm <- lapply(dismats_mrm, function(x) as.dist((x - mean(x)) / sd(x)))
# run mrm
set.seed(64511018)
MRM.results <- list()
for(d in dists){
  MRM.results[[d]] <- MRM(dismats_mrm[[paste0("sal.", d)]] ~ covar_dists_d$host_dist + covar_dists_d$geog_dist + covar_dists_d$clim_dist + covar_dists_d$envm_dist + covar_dists_d$bdlo_dist, nperm = 10000, mrank = TRUE)
}
mrm.summary <- sum.mrm(MRM.results, covar_names = names(covar_dists))
write.csv(mrm.summary$model.tests, file = "results/9.Phylosymbiosis/mrm.model.test.csv", row.names = F)
write.csv(mrm.summary$coef, file = "results/9.Phylosymbiosis/mrm.coef.csv", row.names = F)

######## 4. IS THE SIGNAL OF PHYLOSYMBIOSIS DRIVEN BY HABITAT AND LOCALITY-SPECIFIC MICROBIAL TAXA? ########
## Make dataset with only taxa in all habitat-locality combinations 
all.LocHab <- global_ASVs(physeq = all_f, col = "Loc_Hab")
sal_all_LocHab <- prune_taxa(taxa_names(all.LocHab), sal_r)
sal_all_LocHab_MD <- data.frame(sample_data(sal_all_LocHab))
# get distances
dismats_g <- list()
dists_g <- c("jaccard", "bray", "unifrac", "wunifrac")
for(d in dists_g){
  if(d == "jaccard"){
    dismats_g[[d]] <- phyloseq::distance(sal_all_LocHab, d, binary = TRUE)
  } else {
    dismats_g[[d]] <- phyloseq::distance(sal_all_LocHab, d)
  }
}
# get median microbiome distances
med.sal.dist.per.spp_g <- list() 
for(d in dists_g){
  mymat <- summarise_mat(x = as.matrix(dismats_g[[d]]), metadata = sal_all_LocHab_MD, col = "Species", fun = "median")
  diag(mymat) <- 0
  mymat <- as.dist(mymat)
  med.sal.dist.per.spp_g[[d]] <- mymat
}
rm(list = c("d", "mymat"))
## get microbiome tree
spp.clust_g <- list()
for(d in dists_g){
  tree <- ape::nj(med.sal.dist.per.spp_g[[d]])
  anc <- mrca(tree)["A. jeffersonianum", "N. viridescens"]
  tree <- root(tree, node = anc)
  spp.clust_g[[d]] <- tree
}
## re-test phylosymbiosis
# individual-level mantel test:
set.seed(61511023)
mantel.ind.results_g <- list()
for(d in dists_g){
  mantel.ind.results_g[[d]] <- vegan::mantel(xdis = covar_dists$host_dist, ydis = dismats_g[[d]], permutations = 10000)
}
mantel.ind.summary_g <- sum.mantel(mantel.ind.results_g)
write.csv(mantel.ind.summary_g, file = "results/9.Phylosymbiosis/ind.mantel.gASVs.csv", row.names = F)
# species-level mantel test:
set.seed(63737408)
mantel.spp.results_g <- list()
for(d in dists_g){
  mantel.spp.results_g[[d]] <- vegan::mantel(xdis = covar_dists_spp$host_dist, ydis = med.sal.dist.per.spp_g[[d]], permutations = 10000)
}
mantel.spp.summary_g <- sum.mantel(mantel.spp.results_g)
write.csv(mantel.spp.summary_g, file = "results/9.Phylosymbiosis/spp.mantel.gASVs.csv", row.names = F)
# cospeciation test:
set.seed(97300012)
cospec.results_g <- list()
for(d in dists_g){
  cospec.results_g[[d]] <-  cospeciation(host.phylo, spp.clust_g[[d]], distance = "RF", method = "permutation", assoc = cbind(host.phylo$tip.label, host.phylo$tip.label), nsim = 10000)
}
cospec.summary_g <- sum.cospec(cospec.results_g)
write.csv(cospec.summary_g, file = "results/9.Phylosymbiosis/RF.perm.test.gASVs.csv", row.names = F)
# MRM
# get scaled covar dists for global ASVs 
covar_dists_g <- lapply(covar_dists, function(x) usedist::dist_subset(as.dist(x), rownames(as.matrix(dismats_g$jaccard))))
covar_dists_d_g <- lapply(covar_dists_g, function(x) as.dist((x - mean(x)) / sd(x)))
dismats_g <- lapply(dismats_g, function(x) as.dist((x - mean(x)) / sd(x)))
# run MRM
set.seed(11046443)
MRM.results_g <- list()
for(d in dists_g){
  MRM.results_g[[d]] <- MRM(dismats_g[[d]] ~ covar_dists_d_g$host_dist + covar_dists_d_g$geog_dist + covar_dists_d_g$clim_dist + covar_dists_d_g$envm_dist + covar_dists_d_g$bdlo_dist, nperm = 10000, mrank = TRUE)
}
mrm.summary_g <- sum.mrm(MRM.results_g, covar_names = names(covar_dists_d_g))
write.csv(mrm.summary_g$model.tests, file = "results/9.Phylosymbiosis/mrm.model.test.gASVs.csv", row.names = F)
write.csv(mrm.summary_g$coef, file = "results/9.Phylosymbiosis/mrm.coef.gASVs.csv", row.names = F)
## test phylosymbiosis within habitat subsets (stream and forest only, since pond only contains 2 species)
# prep data
hab_hdi <- list() # individual-level host phylogenetic distance
hab_mdi <- list() # individual-level microbiome distances
for(h in c("Stream", "Forest")){
  # metadata
  md <- sal.metadata[which(sal.metadata$Habitat == h),] 
  # sample list
  ind <- md$sample
  # spp list
  spp <- sort(unique(md$Species)) 
  # host phylo distance matrix
  # species level
  phy <- ape::drop.tip(phy = host.phylo, setdiff(host.phylo$tip.label, spp)) 
  hds <- phy %>% ape::cophenetic.phylo()
  # individual level
  template.dm <- as.matrix(usedist::dist_subset(dismats$sal.ja, ind))
  hdi <- expand_mat(template.dm, hds, metadata = md, col = "Species")
  # per-individual microbiome distance matrices 
  mdi <- list()
  for(d in dists){
    dm <-  dismats[[paste0("sal.",d)]]
    my.mdi <- usedist::dist_subset(dm, ind)
    mdi[[d]] <- my.mdi
  }
  # assign to outputs
  hab_hdi[[h]] <- hdi
  hab_mdi[[h]] <- mdi
}
rm(list = c("d", "dm", "h", "hdi", "hds", "ind", "md", "mdi", "my.mdi", "phy", "spp", "template.dm"))
## re-test phylosymbiosis with individual-level mantel test ():
set.seed(4576434)
mantel.ind.results.hab <- list()
mantel.ind.summary.hab <- list()
for(h in c("Stream", "Forest")){
  mantel.ind.results.hab[[h]] <- list()
  for(d in dists){
    mantel.ind.results.hab[[h]][[d]] <- vegan::mantel(xdis = hab_hdi[[h]], ydis = hab_mdi[[h]][[d]], permutations = 10000)
  }
  mantel.ind.summary.hab[[h]] <- sum.mantel(mantel.ind.results.hab[[h]])
  write.csv(mantel.ind.summary.hab[[h]], file = paste0("results/9.Phylosymbiosis/ind.mantel.", h, ".csv"), row.names = F)
}
rm(list = c( "d", "h"))
# MRM
# get scaled covar dists for habitat subsets
covar_dists.hab <- list()
micro_dists.hab <- list()
for(h in c("Stream", "Forest")){
  cvd <- lapply(covar_dists, function(x) usedist::dist_subset(as.dist(x), rownames(as.matrix(hab_mdi[[h]]$ja))))
  covar_dists.hab[[h]] <- lapply(cvd, function(x) as.dist((x - mean(x)) / sd(x)))
  micro_dists.hab[[h]] <- lapply(hab_mdi[[h]], function(x) as.dist((x - mean(x)) / sd(x)))
}
# run MRM
set.seed(11046443)
MRM.results.hab <- list()
MRM.summary.hab <- list()
for(h in c("Stream", "Forest")){
  MRM.results.hab[[h]] <- list()
  for(d in dists){
    MRM.results.hab[[h]][[d]] <- MRM(micro_dists.hab[[h]][[d]] ~ covar_dists.hab[[h]]$host_dist + covar_dists.hab[[h]]$geog_dist + covar_dists.hab[[h]]$clim_dist + covar_dists.hab[[h]]$envm_dist + covar_dists.hab[[h]]$bdlo_dist, nperm = 10000, mrank = TRUE)
  }
  MRM.summary.hab[[h]] <- sum.mrm(MRM.results.hab[[h]], covar_names = names(covar_dists.hab[[h]]))
  write.csv(MRM.summary.hab[[h]]$model.tests, file = paste0("results/9.Phylosymbiosis/mrm.model.test.", h, ".csv"), row.names = F)
  write.csv(MRM.summary.hab[[h]]$coef, file = paste0("results/9.Phylosymbiosis/mrm.coef.", h, ".csv"), row.names = F)
}
######## 5. DOES PHYLOSYMBIOSIS DRIVE THE DISTRIBUTION OF PATHOGEN-PROTECTIVE MICROBIAL TAXA? ########
## load antiBd
antiBd <- read.csv("results/5.Antifungal_DB/blast.filt.uniq.csv", header = T)
## Make dataset with only taxa in all habitat-locality combinations 
sal_r_abd <- prune_taxa(antiBd$qseqid, sal_r)
sal_r_abd <- prune_samples(sample_sums(sal_r_abd) > 0 , sal_r_abd)
sal_r_abd_MD <- data.frame(sample_data(sal_r_abd))
# get distances
dismats_a <- list()
dists_a <- c("jaccard", "bray", "unifrac", "wunifrac")
for(d in dists_a){
  if(d == "jaccard"){
    dismats_a[[d]] <- phyloseq::distance(sal_r_abd, d, binary = TRUE)
  } else {
    dismats_a[[d]] <- phyloseq::distance(sal_r_abd, d)
  }
}
# get median microbiome distances
med.sal.dist.per.spp_a <- list() 
for(d in dists_a){
  mymat <- summarise_mat(x = as.matrix(dismats_a[[d]]), metadata = sal_r_abd_MD, col = "Species", fun = "median", na.rm=T)
  diag(mymat) <- 0
  mymat <- as.dist(mymat)
  med.sal.dist.per.spp_a[[d]] <- mymat
}
rm(list = c("d", "mymat"))
## get microbiome tree
spp.clust_a <- list()
for(d in dists_a){
  tree <- ape::nj(med.sal.dist.per.spp_a[[d]])
  anc <- mrca(tree)["A. jeffersonianum", "N. viridescens"]
  tree <- root(tree, node = anc)
  spp.clust_a[[d]] <- tree
}
# get host phylogenetic dist per ind
host_dist_abd <- expand_mat(as.matrix(dismats_a$jaccard), as.matrix(covar_dists_spp$host_dist), metadata = sal_r_abd_MD, col = "Species")
# individual-level mantel test
set.seed(61511023)
mantel.ind.results_a <- list()
for(d in dists_a){
  mantel.ind.results_a[[d]] <- vegan::mantel(xdis = host_dist_abd, ydis = dismats_a[[d]], permutations = 10000)
}
mantel.ind.summary_a <- sum.mantel(mantel.ind.results_a)
write.csv(mantel.ind.summary_a, file = "results/9.Phylosymbiosis/ind.mantel.antiBd.csv", row.names = F)
# species-level mantel test:
set.seed(63737408)
mantel.spp.results_a <- list()
for(d in dists_a){
  mantel.spp.results_a[[d]] <- vegan::mantel(xdis = covar_dists_spp$host_dist, ydis = med.sal.dist.per.spp_a[[d]], permutations = 10000)
}
mantel.spp.summary_a <- sum.mantel(mantel.spp.results_a)
write.csv(mantel.spp.summary_g, file = "results/9.Phylosymbiosis/spp.mantel.antiBd.csv", row.names = F)
# cospeciation test:
set.seed(83084)
cospec.results_a <- list()
for(d in dists_a){
  cospec.results_a[[d]] <-  cospeciation(host.phylo, spp.clust_a[[d]], distance = "RF", method = "permutation", assoc = cbind(host.phylo$tip.label, host.phylo$tip.label), nsim = 10000)
}
cospec.summary_a <- sum.cospec(cospec.results_a)
write.csv(cospec.summary_a, file = "results/9.Phylosymbiosis/RF.perm.test.antiBd.csv", row.names = F)
# MRM
# get scaled covar dists for antibd
covar_dists_a <- lapply(covar_dists, function(x) usedist::dist_subset(as.dist(x), rownames(as.matrix(dismats_a$jaccard))))
covar_dists_d_a <- lapply(covar_dists_a, function(x) as.dist((x - mean(x)) / sd(x)))
dismats_a <- lapply(dismats_a, function(x) as.dist((x - mean(x)) / sd(x)))
# run MRM
set.seed(11046443)
MRM.results_a <- list()
for(d in dists_a){
  MRM.results_a[[d]] <- MRM(dismats_a[[d]] ~ covar_dists_d_a$host_dist + covar_dists_d_a$geog_dist + covar_dists_d_a$clim_dist + covar_dists_d_a$envm_dist + covar_dists_d_a$bdlo_dist, nperm = 10000, mrank = TRUE)
}
mrm.summary_a <- sum.mrm(MRM.results_a, covar_names = names(covar_dists_d_a))
write.csv(mrm.summary_a$model.tests, file = "results/9.Phylosymbiosis/mrm.model.test_antiBd.csv", row.names = F)
write.csv(mrm.summary_a$coef, file = "results/9.Phylosymbiosis/mrm.coef_antiBd.csv", row.names = F)
## plot MRM coefficients for all subsets
# set par
par(mar=c(10, 5.1, 2, 2.1), mfrow = c(3,2))
# plot
mrm.sums <- list("All" = mrm.summary,
                 "Global ASVs" = mrm.summary_g,
                 "Stream" = MRM.summary.hab$Stream,
                 "Forest" = MRM.summary.hab$Forest,
                 "Bd-inhibitory ASVs" = mrm.summary_a)
for(n in names(mrm.sums)){
  my.mrm <- mrm.sums[[n]]
  # plot mrm coefficients 
  xpos <- barplot(as.matrix(my.mrm$coef[,c("coeff.host_dist", "coeff.geog_dist", "coeff.clim_dist", "coeff.envm_dist", "coeff.bdlo_dist")]), beside = T, names.arg = rep("", 5), las = 1, ylab = "Regression coeff.", col = grey.colors(4), cex.lab = 1.5, ylim = c(0, 0.7), main = n, cex.main = 2) 
  mtext(at = c(3,8,13,18, 23), side = 1, line = 4, text = c("Host\nphylogeny", "Geographic\ndistance", "Climate\ndistance", "Env. microbiome\ndistance", "Bd\nload"), cex = 1)
  # add stars
  ypos  <- c(my.mrm$coef$coeff.host_dist, my.mrm$coef$coeff.geog_dist,my.mrm$coef$coeff.clim_dist, my.mrm$coef$coeff.envm_dist, my.mrm$coef$coeff.bdlo_dist)
  ypos <- ifelse(ypos > 0, ypos + 0.03, 0.03)
  pvals  <- c(my.mrm$coef$P.host_dist, my.mrm$coef$P.geog_dist, my.mrm$coef$P.clim_dist, my.mrm$coef$P.envm_dist, my.mrm$coef$P.bdlo_dist)
  stars <- stars.pval(pvals)
  text(x = xpos, y = ypos, labels = stars)
}
plot.new()
legend("center", legend = c("Jaccard", "Bray-Curtis", "UniFrac", "Weighted UniFrac"), fill = grey.colors(4), bty = "n", cex = 2, ncol = 1)


######## 6. DOES HOST-MICROBE CO-SPECIATION CONTRIBUTE TO PHYLOSYMBIOSIS? ########
## cluster ASVs
# output ASVs for vsearch & swarm
sw_IDs <- paste(rownames(otu_table(sal_r)),rowSums(otu_table(sal_r)), sep = "_")
vs_IDs <- paste(rownames(otu_table(sal_r)), ";size=", rowSums(otu_table(sal_r)), sep = "")
seqs <- refseq(sal_r)
names(seqs) <- vs_IDs
writeXStringSet(seqs, "results/9.Phylosymbiosis/for_vs.fasta")
names(seqs) <- sw_IDs
writeXStringSet(seqs, "results/9.Phylosymbiosis/for_sw.fasta")
## clustering
# swarm and vsearch should be installed and available in the path
threads <- 6
clust <- list()
min_clust_size <- 3
# swarm clustering
swarmpath <- "bin/swarm"
system(paste(swarmpath, "-t", threads, "results/9.Phylosymbiosis/for_sw.fasta > results/9.Phylosymbiosis/swarms.txt")) # run swarm
clust[["swarm"]] <- readLines("results/9.Phylosymbiosis/swarms.txt") %>% # read output 
  strsplit(split = " ") %>% # split clusters
  list.filter(., length(.) >= min_clust_size ) %>% # remove small clusters
  lapply(X = ., FUN = function(x) gsub(pattern = "_.*", replacement = "", x = x)) # remove abundance from ASV names
# vsearch clustering
vspath <- "bin/vsearch"
for(id in c(0.95, 0.97, 0.99)){
  system(paste(vspath," -cluster_size results/9.Phylosymbiosis/for_vs.fasta --id ", id, " --uc results/9.Phylosymbiosis/clusters_", id, ".uc --threads ", threads, sep = "")) # run vsearch
  clust[[paste0("vs_",id)]] <- read.table(paste0("results/9.Phylosymbiosis/clusters_", id, ".uc")) %>% # read results
    filter(V1 %in% c("S", "H")) %>% # keep only centroids and hits
    arrange(V2) %>% # sort by cluster ID
    group_split(V2) %>% # split by cluster
    lapply(., function(x) pull(x, V9)) %>% # get ASV names
    list.filter(., length(.) >= min_clust_size ) %>% # remove small clusters
    lapply(X = ., FUN = function(x) gsub(pattern = ";.*", replacement = "", x = x)) # remove abundance from ASV names
}
## cospeciation analysis
# merge species
merged_spp_sal_r <- merge_samples(sal_r, group = "Species")
# rarify to ensure even depth
merged_spp_sal_r <- rarefy_even_depth(merged_spp_sal_r, rngseed = 9864210, trimOTUs = T)
# run parafit
set.seed(66839)
parafit_results <- list()
for(cl in names(clust)){
  parafit_results[[cl]] <- run.global_parafit(host_tree = host.phylo, ps = merged_spp_sal_r, clusters = clust[[cl]], min.size = 3, nperm = 10000)
  parafit_results[[cl]]$Padj <- p.adjust(parafit_results[[cl]]$P, method = "fdr")
}
# get average cluster taxonomy
clust.tax <- list()
for(cl in names(clust)){
  clust.tax[[cl]] <- do.call(rbind, sapply(1:length(clust[[cl]]), function(x) consensus_taxonomy(tax_table(sal_r)[clust[[cl]][[x]]], name = paste0("clust.",x)), simplify = F))
}
#save(clust.tax, file = "results/9.Phylosymbiosis/clust.tax.RData")
# add taxonomy and ASV IDs
for(cl in names(clust)){
  parafit_results[[cl]] <- cbind(parafit_results[[cl]], clust.tax[[cl]][,-1])
  parafit_results[[cl]] <- cbind(cluster_method = cl, parafit_results[[cl]], ASVs = unlist(lapply(clust[[cl]], function(x) paste(x, collapse = ";"))))
}
parafit_all <- do.call(rbind, parafit_results)
# save results
write.csv(parafit_all, file = "results/9.Phylosymbiosis/parafit.csv", row.names = F)

