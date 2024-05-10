library(phyloseq)
library(ape)
library(raster) # for getData
library(geodist) # for geodist
library(ade4) # for dudi.pca
source("code/phylosymbiosis_funcs.R") # custom functions

## load phyloseq objects 
load("results/4.Rarefy_and_subset/env_raref.RData")
load("results/4.Rarefy_and_subset/sal_raref.RData")
### load microbiome distance
load("results/7.Beta_diversity/dismats.RData")
# add combined locality-habitat variable
sample_data(sal_r)$Loc_Hab <- gsub(" ", ".", paste(sample_data(sal_r)$Locality, sample_data(sal_r)$Habitat, sep = "_"))
sample_data(sal_r)$sample <- sample_names(sal_r)
sample_data(env_r)$Loc_Hab <- gsub(" ", ".", paste(sample_data(env_r)$Locality, sample_data(env_r)$Habitat, sep = "_"))
# extract metadata
sal.metadata <- data.frame(sample_data(sal_r))
env.metadata <- data.frame(sample_data(env_r))
### skin microbiome distance covariates
covar_dists <- list()
## host phylogenetic distance
host.phylo <- read.tree("data/timetree_31JUL23/host.timetree.31JUL23_short.nwk")
host.phylo$tip.label <- gsub("_", " ", host.phylo$tip.label)
# get distance matrix
host_dist_sp <- cophenetic(host.phylo)
# expand
covar_dists$host_dist <- expand_mat(as.matrix(dismats$sal.ja), as.matrix(host_dist_sp), metadata = sal.metadata, col = "Species")
## geographic distance
coords <- setNames(sal.metadata[,c("longitude", "latitude")], c("x", "y"))
geog_dist <- geodist(coords, measure = "vincenty")
dimnames(geog_dist) <- list(rownames(sal.metadata),rownames(sal.metadata))
covar_dists$geog_dist <- geog_dist
## climatic distance
# get bioclim data
bioclim <- raster::getData("worldclim",var="bio", res = 0.5, lon = sal.metadata[1,c("longitude")], lat = sal.metadata[1,c("latitude")])
# extract values for sampling sites
bc_vals <- extract(bioclim, coords)
rownames(bc_vals) <- rownames(coords)
bc_vals <- as.data.frame(bc_vals)
rm(bioclim)
# add elevation
bc_vals$elev <- sal.metadata$elevation
# run pca
bc_pca <- dudi.pca(bc_vals, scannf = FALSE, nf = 4)
# get climate distance matrix
covar_dists$clim_dist <- as.matrix(dist(bc_pca$li))
## environmental microbiome distance
# first get the mean bray-curtis distance between environmental samples in each locality-habitat combination
envm_dist <- summarise_mat(x = as.matrix(dismats$env.bc), metadata = env.metadata, col = "Loc_Hab", fun = "mean")
# set the diagonal to 0 so there is 0 distance between samples from the same locality and habitat. While samples vary within a habitat-locality, there is no 1-to-1 correspondence between environmental and salamander samples so there would otherwise just be the same average between all members of the same group.
diag(envm_dist) <- 0
# then expand the matrix to give mean environmental microbiome distance between each pair of salamander samples
envm_dist <- expand_mat(as.matrix(dismats$sal.ja), envm_dist, metadata = sal.metadata, col = "Loc_Hab")
covar_dists$envm_dist <- envm_dist 
## Bd load distance
bdload <- sal.metadata$Bd_load
bdload <- ifelse(is.na(bdload),0,bdload)
bdload <- log10(bdload+1)
names(bdload) <- rownames(sal.metadata)
bdlo_dist <- as.matrix(dist(bdload)) 
covar_dists$bdlo_dist <- bdlo_dist
# check all covariates are in the same order
for(n in names(covar_dists)) {print(all(rownames(covar_dists[[n]]) == rownames(covar_dists$host_dist))) ; print(all(colnames(covar_dists[[n]]) == colnames(covar_dists$host_dist)))}
# get covariates in same order
covar_dists_spp <- lapply(covar_dists, FUN = function(x) as.dist(summarise_mat(x, metadata = sal.metadata, col = "Species", fun = "mean", diagzero = TRUE)))
# save
dir.create("results/8.Environmental_distances/")
save(list = c("covar_dists","covar_dists_spp"), file = "results/8.Environmental_distances/env.dists.RData")
