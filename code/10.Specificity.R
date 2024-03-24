# load pkgs
library(specificity)
library(metacoder)
library(phyloseq)
library(dplyr)
library(colorspace)
library(gtools)
source("code/phylosymbiosis_funcs.R") # custom functions
dir.create("results/10.Specificity/")
# load data
load("results/4.Rarefy_and_subset/sal_raref.RData")
## load env distances 
load("results/8.Environmental_distances/env.dists.RData")
# get OTU tab
sal.otutable <-  t(as(otu_table(sal_r), "matrix"))
# transform to proportional abundance
#sal.otutableP <- prop_abund(sal.otutable)
# remove taxa in less than 10 samples
sal.otutable_ovr10 <- occ_threshold(sal.otutable, threshold=10)
# set cores
spec_ncores <- 6
# run
spec_all <- list()
set.seed(44532)
for(e in names(covar_dists)){
  spec_all[[e]] <- phy_or_env_spec(sal.otutable_ovr10, env = covar_dists[[e]], n_sim=10000, n_cores = spec_ncores, tails = 2)
}
spec_comb <- do.call(cbind, spec_all)
spec_comb$AntiBd <- names(spec_comb) %in% AntiBd$qseqid
# plot correlations between specificities
names(spec_all) <- c("Host Phylogeny", "Geographic Distance", "Climate Distance", "Environmental Microbiome Distance")
plot_pairwise_spec(spec_all, cor_red_lim = 1, label_cex = 1.2, point_cex = 0.5)
# save
save(list = c("spec_all", "spec_comb"), file = "results/10.Specificity/specificity_res.RData")
## plot number of significant specificity results beneath MRM results
# load mrm
mrm_coef <- read.csv("results/9.Phylosymbiosis/mrm.coef.csv")
# save par
OLDPAR <- par()
# set par
par(mfrow = c(2,1))
par(mar=c(1, 5.1, 1, 2.1))
# plot mrm
bp <- barplot(as.matrix(mrm_coef[,c("coeff.host_dist", "coeff.geog_dist", "coeff.clim_dist", "coeff.envm_dist")]), beside = T, names.arg = rep("", 4), las = 2, ylab = "Regression coeff.", col = grey.colors(4), cex.lab = 1.5, ylim = c(0, 0.7)) #c("Host phylogeny", "Geographic distance", "Climate distance", "environmental microbiome distance")
legend("topright", legend = c("Jaccard", "Bray-Curtis", "UniFrac", "Weighted UniFrac"), fill = grey.colors(4), bty = "n", cex = 1.5, ncol = 2)
# add stars
xpos <- bp
ypos  <- c(mrm_coef$coeff.host_dist, mrm_coef$coeff.geog_dist, mrm_coef$coeff.clim_dist, mrm_coef$coeff.envm_dist)
ypos <- ifelse(ypos > 0, ypos + 0.03, 0.03)
pvals  <- c(mrm_coef$P.host_dist, mrm_coef$P.geog_dist, mrm_coef$P.clim_dist, mrm_coef$P.envm_dist)
stars <- stars.pval(pvals)
text(x = xpos, y = ypos, labels = stars)
# plot specificity
Nspec <- unlist(lapply(spec_all, function(x) length(which(x$Pval <= 0.05 & x$Spec < 0))))
Ncosm <- unlist(lapply(spec_all, function(x) length(which(x$Pval <= 0.05 & x$Spec > 0))))
barplot(rbind(Nspec, Ncosm), col =c("gray", "black"), border = NA, ylab = "N significant ASVs", cex.lab = 1.5)
legend("top", legend = c("Specific", "Cosmopolitan"), fill =c("gray", "black"), bty = "n", cex = 1.5)
# restore default par
par(OLDPAR)
## plot heat tree
# get taxonomy and reformat
tax <- as(tax_table(sal_r), "matrix") 
lin <- paste0("r__Root;p__", tax[,"Phylum"],
              ";c__", tax[,"Class"],
              ";o__", tax[,"Order"],
              ";f__", tax[,"Family"],
              ";g__", tax[,"Genus"],
              ";s__", tax[,"Species"],
              ";a__", tax[,"ASV"])
otuid <- tax[,"ASV"]
# make input df
tt <- data.frame(otu_id = otuid, lineage = lin)
# add specificity analysis
tt <- tt[rownames(spec_comb),]
tt <- cbind(tt, spec_comb)
# get otu_tab
otutab <- data.frame(otu_table(sal_r))
otutab <- otutab[rownames(spec_comb),]
# combine
tt <- cbind(tt, otutab)
# make tax_data object
x = parse_tax_data(tt, class_cols = "lineage", class_sep = ";",
                   class_key = c(tax_rank = "taxon_rank", tax_name = "taxon_name"),
                   class_regex = "^(.+)__(.+)$")
x$data$host_spec <- unlist(lapply(obs(x, data = 'tax_data'), function(i)  mean(x$data$tax_data$host_dist.Spec[i])))
x$data$envm_spec <- unlist(lapply(obs(x, data = 'tax_data'), function(i)  mean(x$data$tax_data$envm_dist.Spec[i])))
x$data$clim_spec <- unlist(lapply(obs(x, data = 'tax_data'), function(i)  mean(x$data$tax_data$clim_dist.Spec[i])))
x$data$geog_spec <- unlist(lapply(obs(x, data = 'tax_data'), function(i)  mean(x$data$tax_data$geog_dist.Spec[i])))


# taxon abundance
x$data$tax_table <- calc_taxon_abund(x, data = "tax_data", cols = sample_names(sal_r), groups = rep("Abundance", nsamples(sal_r)))
# palette
coolwarm_hcl <- rev(diverging_hcl(11, h = c(250, 10), c = 100, l = c(37, 88), power = c(0.7, 1.7)))
# plot host phylogenetic specificity (only order level and above)
set.seed(1234)
x %>% 
  metacoder::filter_taxa(taxon_ranks == "o", supertaxa = T) %>% 
  heat_tree(node_label = ifelse(!(taxon_ranks %in% c("p", "c", "o")) | taxon_names == "NA", "", taxon_names), 
            node_size = Abundance, 
            node_color = host_spec, 
            node_color_range = coolwarm_hcl, 
            layout = "davidson-harel", 
            initial_layout = "reingold-tilford", 
            node_color_interval = c(-1, 1),
            edge_color_interval = c(-1, 1), 
            node_color_trans = "linear",
            node_size_trans = "linear",
            node_label_size_trans ="area",
            node_color_axis_label = "Specificity",
            output_file = "results/10.Specificity/host.spec.heattree.o.pdf")
# plot env microbiome specificity
set.seed(1234)
x %>% 
  metacoder::filter_taxa(taxon_ranks == "o", supertaxa = T) %>% 
  heat_tree(node_label = ifelse(!(taxon_ranks %in% c("p", "c", "o")) | taxon_names == "NA", "", taxon_names), 
            node_size = Abundance, 
            node_color = envm_spec, 
            node_color_range = coolwarm_hcl, 
            layout = "davidson-harel", 
            initial_layout = "reingold-tilford", 
            node_color_interval = c(-1, 1),
            edge_color_interval = c(-1, 1),
            node_color_trans = "linear",
            node_size_trans = "linear",
            node_label_size_trans ="area",
            node_color_axis_label = "Specificity",
            output_file = "results/10.Specificity/envm.spec.heattree.o.pdf")
# plot geographic specificity
set.seed(1234)
x %>% 
  metacoder::filter_taxa(taxon_ranks == "o", supertaxa = T) %>% 
  heat_tree(node_label = ifelse(!(taxon_ranks %in% c("p", "c", "o")) | taxon_names == "NA", "", taxon_names), 
            node_size = Abundance, 
            node_color = geog_spec, 
            node_color_range = coolwarm_hcl, 
            layout = "davidson-harel", 
            initial_layout = "reingold-tilford", 
            node_color_interval = c(-1, 1),
            edge_color_interval = c(-1, 1),
            node_color_trans = "linear",
            node_size_trans = "linear",
            node_label_size_trans ="area",
            node_color_axis_label = "Specificity",
            output_file = "results/10.Specificity/geog.spec.heattree.o.pdf")
# plot climate specificity
set.seed(1234)
x %>% 
  metacoder::filter_taxa(taxon_ranks == "o", supertaxa = T) %>% 
  heat_tree(node_label = ifelse(!(taxon_ranks %in% c("p", "c", "o")) | taxon_names == "NA", "", taxon_names), 
            node_size = Abundance, 
            node_color = clim_spec, 
            node_color_range = coolwarm_hcl, 
            layout = "davidson-harel", 
            initial_layout = "reingold-tilford", 
            node_color_interval = c(-1, 1),
            edge_color_interval = c(-1, 1),
            node_color_trans = "linear",
            node_size_trans = "linear",
            node_label_size_trans ="area",
            node_color_axis_label = "Specificity",
            output_file = "results/10.Specificity/clim.spec.heattree.o.pdf")
## plot envm vs host spec, with known anti-Bd highlighted
# combine Specificity results and taxonomy
tt <- tax_table(sal_r)
tt <- data.frame(tt[rownames(spec_comb),])
df <- cbind(tt,spec_comb)
rm(tt)
# get antiBD 
AntiBd <- read.csv("results/5.Antifungal_DB/blast.filt.uniq.csv", header = T)
df$AntiBd <- df$ASV %in% AntiBd$qseqid
# save
write.csv(df, file = "results/10.Specificity/Specificity_results.csv", row.names = F)
# combine rare clades into "other"
df.mod <- combn_rare_class(df, min.p = 20, min.c = 10)
df.mod.abd <- df.mod[which(df.mod$AntiBd == T),]
# colours
phy.cols <- c("forestgreen", "goldenrod", "cornflowerblue", "lightcoral", "purple", "black")
names(phy.cols) <- c("Acidobacteria", "Actinobacteria", "Bacteroidetes",  "Proteobacteria", "Verrucomicrobia", "Other")
phy.cols.ns <-  adjust_transparency(phy.cols, 0.3)
names(phy.cols.ns) <- names(phy.cols)
# par
par(mar = c(5.1, 6.1, 4.1, 2.1))
# plot
plot(df.mod$host_dist.Spec, df.mod$envm_dist.Spec, pch = 19, col = phy.cols.ns[df.mod$Phylum], cex = 2, xlab = "Host phylogeny Spec.", ylab = "Environmental microbiome Spec.",cex.axis = 1.5, cex.lab = 2, xlim = c(-1, 1), ylim = c(-1, 1), bty = "n")
abline(v=0)
abline(h=0)
points(df.mod.abd$host_dist.Spec, df.mod.abd$envm_dist.Spec, bg = phy.cols[df.mod.abd$Phylum], cex = 3, pch = 23, lwd = 1, col = "black")
legend("bottomright", legend = names(phy.cols), fill = phy.cols, cex = 1.5, bty = "n")
## how many species are low spec anti-Bd taxa in?
least_spec_abd_otus <- df.mod.abd[which(df.mod.abd$host_dist.Spec > 0 & df.mod.abd$envm_dist.Spec > 0),]
# by species
merged_spp_sal_r <- merge_samples(sal_r, group = "Species")
otu_spp <- data.frame(t(otu_table(merged_spp_sal_r)))
otu_spp <- otu_spp[least_spec_abd_otus$ASV,]
ifelse(otu_spp > 0, T, F)
# by habitat
merged_hab_sal_r <- merge_samples(sal_r, group = "Habitat")
otu_hab <- data.frame(t(otu_table(merged_hab_sal_r)))
otu_hab <- otu_hab[least_spec_abd_otus$ASV,]
ifelse(otu_hab > 0, T, F)
