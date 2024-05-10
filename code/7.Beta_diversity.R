library(phyloseq)
library(ggplot2)
library(ape)
library(gridExtra)
library(GUniFrac)
library(usedist)
library(gtools)
source("code/phylosymbiosis_funcs.R") # custom functions
# load data
load("results/4.Rarefy_and_subset/all_raref.RData")
load("results/4.Rarefy_and_subset/sal_raref.RData")
load("results/4.Rarefy_and_subset/env_raref.RData")
metadata <- data.frame(sample_data(all_r))

dists <- c("jaccard", "bray", "unifrac", "wunifrac")
dists.s <- c("ja", "bc", "uu", "wu") ; names(dists.s) <- dists

dismats <- list()
DS <- c("all", "env", "sal")
for(D in DS){
  for(d in dists){
    cat(d, D, "\n")
    if(d == "jaccard"){
      dismats[[paste(D, dists.s[[d]], sep= ".")]] <-
        distance(get(paste(D, "r", sep = "_")), d, binary = TRUE)
    } else {
      dismats[[paste(D, dists.s[[d]], sep= ".")]] <-
        distance(get(paste(D, "r", sep = "_")), d)
    }
  }
}
dir.create("results/7.Beta_diversity")
save(dismats, file = "results/7.Beta_diversity/dismats.RData")
#NMDS
set.seed(290723)
nmds <- list()
for(D in DS){
  for(d in dists){
    cat(d, D, "\n")
    # get inputs
    ps <- get(paste(D, "r", sep = "_"))
    dm <- dismats[[paste(D, dists.s[[d]], sep = ".")]]
    # ordination
    ord <- ordinate(ps,
                    method = "NMDS",
                    distance = dm)
    # add to output
    nmds[[paste(D, dists.s[[d]], sep= ".")]] <- ord
  }
}
# colours
sppcol <- c("black","olivedrab","forestgreen","#007054","lightblue","gray","gold","#F2A181","#FFB0C9","darkmagenta")
names(sppcol) <- sort(unique(metadata$Species))
habcol <- c("#009E73", "#0072B2", "#56B4E9")
names(habcol) <- sort(unique(metadata$Habitat))
loccol <- c("#E69F00", "#D55E00", "#CC79A7")
names(loccol) <- sort(unique(metadata$Locality))
# ordination plots
ord.plots <- list()
D <- "all"
for(d in dists){
  cat(d, D, "\n")
  # get inputs
  ps <- get(paste(D, "r", sep = "_"))
  ord <- nmds[[paste(D, dists.s[[d]], sep= ".")]]
  scaleFUN <- function(x) sprintf("%.2f", x)
  # habitat plot
  hab.plot <- plot_ordination(all_r, ord, shape = "Sample_Type", color = "Habitat") +  
    scale_shape_manual(values = c(4, 19)) + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "white", color = "black")) +
    scale_y_continuous(labels=scaleFUN) +
    scale_x_continuous(labels=scaleFUN) +
    scale_color_manual(values = habcol) + 
    theme(aspect.ratio=1)
  # species plot
  spp.plot <- plot_ordination(all_r, ord, shape = "Sample_Type", color = "Species") +  
    scale_shape_manual(values = c(4, 19)) + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(),
          axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "white", color = "black")) +
    scale_y_continuous(labels=scaleFUN) +
    scale_x_continuous(labels=scaleFUN) +
    scale_color_manual(values = sppcol) + 
    theme(aspect.ratio=1)
  # locality plot
  loc.plot <- plot_ordination(all_r, ord, shape = "Sample_Type", color = "Locality") +  
    scale_shape_manual(values = c(4, 19)) + 
    theme(legend.position = "none", axis.title.x = element_blank(), axis.title.y = element_blank(), 
          axis.text.x=element_text(colour="black"), axis.text.y=element_text(colour="black"),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          panel.background = element_rect(fill = "white", color = "black")) +
    scale_y_continuous(labels=scaleFUN) +
    scale_x_continuous(labels=scaleFUN) +
    scale_color_manual(values = loccol) + 
    theme(aspect.ratio=1)
  # add to output
  ord.plots[[paste(D, dists.s[[d]], "hab", sep= ".")]] <- hab.plot
  ord.plots[[paste(D, dists.s[[d]], "spp", sep= ".")]] <- spp.plot
  ord.plots[[paste(D, dists.s[[d]], "loc", sep= ".")]] <- loc.plot
}

# MAIN FIGURE
grid.arrange(ord.plots$all.bc.spp + ggtitle("A Host species") + theme(plot.title = element_text(size = 20)),
          ord.plots$all.bc.hab + ggtitle("B Habitat") + theme(plot.title = element_text(size = 20)),
          ord.plots$all.bc.loc + ggtitle("C Locality") + theme(plot.title = element_text(size = 20)),
          layout_matrix = matrix(1:3, ncol=1), respect = T) 
# make ones with legend to cut + paste on
# function to make italic legend for species
make.italic <- function(x) as.expression(lapply(x, function(y) bquote(italic(.(y)))))
grid.arrange(ord.plots$all.bc.spp + ggtitle("(a) Host Species") + theme(plot.title = element_text(size = 20), legend.position = "right", legend.key = element_blank(), legend.text.align = 0) + scale_color_manual(values = sppcol, limits = names(sppcol), labels = make.italic(names(sppcol))) + labs(shape = "Sample Type"),
             ord.plots$all.bc.hab + ggtitle("(b) Habitat") + theme(plot.title = element_text(size = 20), legend.position = "right", legend.key = element_blank()),
             ord.plots$all.bc.loc + ggtitle("(c) Locality") + theme(plot.title = element_text(size = 20), legend.position = "right", legend.key = element_blank()),
             layout_matrix = matrix(1:3, ncol=1)) 
# SUPPLEMENTARY FIGURE
grid.arrange(ord.plots$all.ja.spp,
             ord.plots$all.ja.hab,
             ord.plots$all.ja.loc,
             ord.plots$all.bc.spp,
             ord.plots$all.bc.hab,
             ord.plots$all.bc.loc,
             ord.plots$all.uu.spp,
             ord.plots$all.uu.hab,
             ord.plots$all.uu.loc,
             ord.plots$all.wu.spp,
             ord.plots$all.wu.hab,
             ord.plots$all.wu.loc,
             layout_matrix = matrix(1:12, ncol=3, byrow = T), respect = T)

# env 
ppG <- list()
# set seed
set.seed(23820801)
PERM <- 10000
# subset data
env.dists <- list(dismats$env.ja, dismats$env.bc, dismats$env.uu, dismats$env.wu)
env.groups <- data.frame(sample_data(env_r))[,c("Locality", "Habitat")]
# general permanova
perm.env <- PermanovaG2(env.dists ~ Locality * Habitat, data = env.groups, permutations = PERM)
perm.env.tab <- pg_output(perm.env)
perm.env.tab
# pairwise permanova
ppG$env.ppG.hab <- pairwise.permG(my.groups = env.groups, my.col = "Habitat", my.form = "~ Locality * Habitat", my.dists = env.dists, perm = PERM)
ppG$env.ppG.loc <- pairwise.permG(my.groups = env.groups, my.col = "Locality", my.form = "~ Locality * Habitat", my.dists = env.dists, perm = PERM)
# sals
sal.dists <- list(dismats$sal.ja, dismats$sal.bc, dismats$sal.uu, dismats$sal.wu)
sal.groups <- data.frame(sample_data(sal_r))[,c("Species", "Locality", "Habitat")]
# forest
# subset data
sal.for.groups <- sal.groups[which(sal.groups$Habitat == "Forest"), c("Species", "Locality")]
sal.for.dists <- sal.dists
for(j in 1:4){
  sal.for.dists[[j]] <- usedist::dist_subset(sal.for.dists[[j]], rownames(sal.for.groups))
}
# general permanova
perm.sal.for <- PermanovaG2(sal.for.dists ~ Locality * Species, data = sal.for.groups, permutations = PERM)
perm.sal.for.tab <- pg_output(perm.sal.for)
perm.sal.for.tab
# pairwise permanova
ppG$sal.for.ppG.spp <- pairwise.permG(my.groups = sal.for.groups, my.col = "Species", my.form = "~ Locality * Species", my.dists = sal.for.dists, perm = PERM)
ppG$sal.for.ppG.loc <- pairwise.permG(my.groups = sal.for.groups, my.col = "Locality", my.form = "~ Locality * Species", my.dists = sal.for.dists, perm = PERM)

# pond
# subset data
sal.pon.groups <- sal.groups[which(sal.groups$Habitat == "Pond"), c("Species", "Locality")]
sal.pon.dists <- sal.dists
for(j in 1:4){
  sal.pon.dists[[j]] <- usedist::dist_subset(sal.pon.dists[[j]], rownames(sal.pon.groups))
}
# general permanova
perm.sal.pon <- PermanovaG2(sal.pon.dists ~ Locality * Species, data = sal.pon.groups, permutations = PERM)
perm.sal.pon.tab <- pg_output(perm.sal.pon)
perm.sal.pon.tab
# pairwise permanova
ppG$sal.pon.ppG.spp <- pairwise.permG(my.groups = sal.pon.groups, my.col = "Species", my.form = "~ Locality * Species", my.dists = sal.pon.dists, perm = PERM)
ppG$sal.pon.ppG.loc <- pairwise.permG(my.groups = sal.pon.groups, my.col = "Locality", my.form = "~ Locality * Species", my.dists = sal.pon.dists, perm = PERM)

# stream
# subset data
sal.str.groups <- sal.groups[which(sal.groups$Habitat == "Stream"), c("Species", "Locality")]
sal.str.dists <- sal.dists
for(j in 1:4){
  sal.str.dists[[j]] <- usedist::dist_subset(sal.str.dists[[j]], rownames(sal.str.groups))
}
# general permanova
perm.sal.str <- PermanovaG2(sal.str.dists ~ Locality * Species, data = sal.str.groups, permutations = PERM)
perm.sal.str.tab <- pg_output(perm.sal.str)
perm.sal.str.tab
# pairwise permanova
ppG$sal.str.ppG.spp <- pairwise.permG(my.groups = sal.str.groups, my.col = "Species", my.form = "~ Locality * Species", my.dists = sal.str.dists, perm = PERM)
ppG$sal.str.ppG.loc <- pairwise.permG(my.groups = sal.str.groups, my.col = "Locality", my.form = "~ Locality * Species", my.dists = sal.str.dists, perm = PERM)

#Nvir
# subset data
sal.nvi.groups <- sal.groups[which(sal.groups$Species == "N. viridescens"), c("Habitat", "Locality")]
sal.nvi.dists <- sal.dists
for(j in 1:4){
  sal.nvi.dists[[j]] <- usedist::dist_subset(sal.nvi.dists[[j]], rownames(sal.nvi.groups))
}
# general permanova
perm.sal.nvi <- PermanovaG2(sal.nvi.dists ~ Locality * Habitat, data = sal.nvi.groups, permutations = PERM)
perm.sal.nvi.tab <- pg_output(perm.sal.nvi)
perm.sal.nvi.tab
# pairwise permanova
ppG$sal.nvi.ppG.hab <- pairwise.permG(my.groups = sal.nvi.groups, my.col = "Habitat", my.form = "~ Locality * Habitat", my.dists = sal.nvi.dists, perm = PERM)
ppG$sal.nvi.ppG.loc <- pairwise.permG(my.groups = sal.nvi.groups, my.col = "Locality", my.form = "~ Locality * Habitat", my.dists = sal.nvi.dists, perm = PERM)

# combine
ppG <- do.call(rbind, ppG)
# correct pvals
ppG$adj.omni.p <- p.adjust(ppG$omni.p, method = "fdr")
ppG$test <- row.names(ppG)
ppG$sig <- stars.pval(ppG$adj.omni.p)
# save
write.csv(perm.env.tab, file = "results/7.Beta_diversity/permanovaG.env.csv")
write.csv(perm.sal.for.tab, file = "results/7.Beta_diversity/permanovaG.sal.for.csv")
write.csv(perm.sal.pon.tab, file = "results/7.Beta_diversity/permanovaG.sal.pon.csv")
write.csv(perm.sal.str.tab, file = "results/7.Beta_diversity/permanovaG.sal.str.csv")
write.csv(perm.sal.nvi.tab, file = "results/7.Beta_diversity/permanovaG.nvi.pon.csv")
write.csv(ppG, file = "results/7.Beta_diversity/pairwise.permanovaG.csv", row.names = F)

