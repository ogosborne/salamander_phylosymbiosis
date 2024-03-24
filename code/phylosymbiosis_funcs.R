##### GENERAL FUNCTIONS FOR MANIPULATING DISTANCE/DISSIMILARITY MATRICES

expand_mat <- function(x, y, metadata, col){
# this function takes a dissimilarity matrix between samples (x), and another dissimilarity matrix between groups made up of those samples (y) and outputs a matrix with the dimensions of (x) containing the values taken from (y)
# x: dissimilarity matrix by sample 
# y: dissimilarity matrix by group 
# metadata: data frame containing one column called 'Sample' (containing identical values to row/column names of x) and another column with the sample groupings (containing identical values to row/column names of y). Additional columns are ignored.
# col: column of metadata which contains comparisons in matrix y
if(any(row.names(x) != colnames(x)) | any(row.names(y) != colnames(y))){
  stop("column and row names must be identical")
}
# make named vector to convert sample to group names
samp2group <- metadata[,col]
names(samp2group) <- rownames(metadata)
# make output matrix
out <- matrix(NA, nrow = nrow(x), ncol = ncol(x), dimnames = list(rownames(x), colnames(x)))
# fill output matrix
for(i in rownames(x)){
  for(j in colnames(x)){
    if(samp2group[i] %in% rownames(y) & samp2group[j] %in% rownames(y)){
      out[i, j] <-  y[samp2group[i], samp2group[j]]
    } else {
      out[i, j] <- NA
    }
  }
}
out
}

summarise_mat <- function(x, metadata, col, fun = "mean", diagzero = FALSE, ...){
  # this function takes a dissimilarity matrix between samples, and summarises samples by a column (col) in the data frame "metadata". By default it returns the mean of all samples for each value, but can use any function which takes a numeric vector and returns a single number.
  if(any(row.names(x) != colnames(x))){
    stop("column and row names must be identical")
  }
  # make named vector to convert sample to group names
  # 
  out <- matrix(NA, nrow = length(unique(metadata[,col])), ncol = length(unique(metadata[,col])), dimnames = list(sort(unique(metadata[,col])), sort(unique(metadata[,col]))))
  # 
  for (i in sort(unique(metadata[,col]))){
    for (j in sort(unique(metadata[,col]))){
      isamps <- rownames(metadata[which(metadata[,col] == i),])
      jsamps <- rownames(metadata[which(metadata[,col] == j),])
      vals <- c(x[isamps, jsamps])
      res <- try(get(fun)(vals, ...))
      if(class(res) == "try-error") res <- NA
      out[i, j] <- res
    }
  }
  if(diagzero == TRUE){
    diag(out) <- 0
  }
  out
}

##### PermanovaG TESTS
# function to process output 
pg_output <- function(pg.res){
  ja <- as.data.frame(pg.res$aov.tab.list[[1]])
  bc <- as.data.frame(pg.res$aov.tab.list[[2]])
  uu <- as.data.frame(pg.res$aov.tab.list[[3]])
  wu <- as.data.frame(pg.res$aov.tab.list[[4]])
  df <- data.frame(
    var = rownames(ja)[1:(nrow(ja)-1)],
    DF = ja$Df[1:(nrow(ja)-1)],
    R2_ja = round(ja$R2[1:(nrow(ja)-1)],2),
    R2_bc = round(bc$R2[1:(nrow(ja)-1)],2),
    R2_uu = round(uu$R2[1:(nrow(ja)-1)],2),
    R2_wu = round(wu$R2[1:(nrow(ja)-1)],2),
    P = c(pg.res$p.tab$omni.p.value, NA)
  )
  df
}
# pairwise PermanovaG function
pairwise.permG <- function(my.groups, my.col, my.form, my.dists, perm = 99999){
  vals <- sort(unique(my.groups[,my.col]))
  combs <- combn(vals, m = 2)
  out.df <- data.frame(matrix(NA, ncol = 8, nrow = ncol(combs), dimnames = list(NULL, c("grp1", "grp2", "DF", "R2.ja", "R2.bc", "R2.uu", "R2.wu", "omni.p"))))
  for(i in 1:ncol(combs)){
    set <- rownames(my.groups[which(my.groups[,my.col] %in% combs[,i]),])
    set.grp <- my.groups[set,]
    my.set.dists <- my.dists
    for(j in 1:4){
      my.set.dists[[j]] <- usedist::dist_subset(my.dists[[j]], set)
    }
    my.other.col <- colnames(my.groups)[which(colnames(my.groups) != my.col)]
    if(length(unique(set.grp[,my.other.col])) < 2){
      out <- PermanovaG2(formula(paste("my.set.dists ~", my.col)), data = set.grp, permutations = perm)
    } else {
      out <- PermanovaG2(formula(paste("my.set.dists", my.form)), data = set.grp, permutations = perm)
    }
    # extract output
    out.df[i,1:2] <- combs[,i]
    out.df[i,3] <- out$aov.tab.list[[1]][my.col,"Df"]
    out.df[i,4:8] <- c(out$aov.tab.list[[1]][my.col,"R2"], out$aov.tab.list[[2]][my.col,"R2"], out$aov.tab.list[[3]][my.col,"R2"], out$aov.tab.list[[4]][my.col,"R2"], out$p.tab[my.col,"omni.p.value"])
  }
  out.df
}
### FUNCTIONS FOR SUMMARISING TEST RESULTS
# Put list of mantel results into a data frame
sum.mantel <- function(x){
  n <- names(x)
  out <- data.frame(matrix(NA, ncol = 3, nrow = length(n), dimnames = list(n, c("name", "R", "P"))))
  out$name <- n
  for(N in 1:length(n)){
    out[N,2:3] <- c(x[[n[[N]]]]$statistic, x[[n[[N]]]]$signif)
  }
  out
} 
# Put list of cospeciation results into a data frame
sum.cospec <- function(x){
  n <- names(x)
  out <- data.frame(matrix(NA, ncol = 5, nrow = length(n), dimnames = list(n, c("name", "obs.d", "mean.null.d", "sd.null.d", "P"))))
  out$name <- n
  for(N in 1:length(n)){
    out[N,2:5] <- c(x[[n[[N]]]]$d, mean(x[[n[[N]]]]$d.null), sd(x[[n[[N]]]]$d.null), x[[n[[N]]]]$P.val)
  }
  out
} 
# Put list of cospeciation results into a data frame
sum.mrm <- function(x, covar_names){
  dists <- names(x)
  ndist <- length(dists)
  nvars <- nrow(x[[1]]$coef)
  outRF <- data.frame(matrix(NA, ncol = 5, nrow = ndist, dimnames = list(dists, c("dist", "R2", "Pval.R2",  "F", "Pval.F"))))
  outCo <- data.frame(matrix(NA, ncol = nvars * 2 + 1, nrow = ndist, dimnames = list(dists, c("dist", "coeff.int", paste0("coeff.", covar_names), "P.int", paste0("P.", covar_names)))))
  outCo$dist <- outRF$dist <- dists
  for(N in 1:ndist){
    outRF[N,2:5] <- c( x[[dists[[N]]]]$r.squared, x[[dists[[N]]]]$F.test)
    outCo[N,2:11] <- c(round(x[[dists[[N]]]]$coef[,1],2), x[[dists[[N]]]]$coef[,2])
  }
  list(model.tests = outRF, coef = outCo)
} 

###### FUNCTION TO PLOT COPHYLO

plot.col.cophylo <- function(t1, t2, col, link.alpha = 1, point.cex = 2.5, link.lwd = 4, link.lty = 1, link.type = "curved"){
  # shorten tip labels to numbers
  tl_n <- 1:length(t1$tip.label)
  names(tl_n) <- t1$tip.label
  t1$tip.label <- tl_n[as.character(t1$tip.label)]
  t2$tip.label <- tl_n[as.character(t2$tip.label)]
  # change palette labels by new tip labels
  names(col) <- tl_n[as.character(names(col))]
  link.col <- phytools::make.transparent(col, link.alpha)
  # make cophylo
  cp <- phytools::cophylo(tr1 = t1, tr2 = t2, rotate = T, methods = c("pre","post"))
  # plot cophylo
  plot(cp, fsize=1E-20, link.lty = link.lty , link.lwd = link.lwd, link.type = link.type, link.col = link.col[as.character(cp$assoc[,1])], edge.lengths = FALSE)
  # add tip markers 
  tiplabels.cophylo(pch = 21, cex = point.cex, bg = col[as.character(cp$trees[[1]]$tip.label)], which = "left")
  tiplabels.cophylo(pch = 21, cex = point.cex, bg = col[as.character(cp$trees[[2]]$tip.label)], which = "right")
}

##### FUNCTION TO GET ONLY ASVS IN ALL ENVS
# get ASVs in all envs
global_ASVs <- function(physeq, col){
  metadata <- data.frame(phyloseq::sample_data(physeq))
  sites <- sort(unique( metadata[,col]))
  all.tax <- list()
  for( s in sites){
    my.ps <- phyloseq::prune_samples(data.frame(phyloseq::sample_data(physeq))[,col] == s, physeq)
    my.ps <- phyloseq::filter_taxa(my.ps, function(x) sum(x) != 0, TRUE)
    all.tax[[s]] <- taxa_names(my.ps)
  }
  global <- Reduce(intersect, all.tax)
  out <- prune_taxa(global, physeq)
  out
}
##### FUNCTION USED IN THE PARAFIT ANALYSIS
consensus_taxonomy <- function(tax.tab, levels = c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), name = "OTU1"){
  tax.tab <- data.frame(tax.tab)[,levels]
  nocons <- which(apply(tax.tab, 2, function(x) length(unique(x))) > 1)
  tax.tab <- tax.tab[1,]
  tax.tab[1,nocons] <- NA
  rownames(tax.tab) <- name
  tax.tab <- cbind(clust = name, tax.tab)
  tax.tab
}

run.global_parafit = function(host_tree, ps, clusters, min.size = 3, nperm = 10000){
  # initiate output
  res.out <- data.frame(matrix(NA, ncol = 3, nrow = length(clusters), dimnames = list(NULL, c("clust", "stat", "P")))) 
  res.out$clust <- paste0("clust.",1:length(clusters))
  #ps.out <- list()
  # progress bar
  pb <- utils::txtProgressBar(max = length(clusters), style = 3)
  # for each cluster
  for(cluster in 1:length(clusters)){
    utils::setTxtProgressBar(pb, cluster)
    # subset phyloseq object
    my.asvs <- clusters[[cluster]]
    my.ps <- try(phyloseq::prune_taxa(phyloseq::tax_table(ps)[,"ASV"] %in% my.asvs, ps), silent = T)
    my.ps <- try(phyloseq::prune_samples(phyloseq::sample_sums(my.ps) > 0, my.ps), silent = T)
    # check subsetted PS still has enough samples
    if(class(my.ps) != "try-error"){
      if(phyloseq::ntaxa(my.ps) >= min.size & phyloseq::nsamples(my.ps) >= min.size){
        # prep input
        my.HP <- phyloseq::otu_table(my.ps) %>% as(., "matrix") %>% t() %>% apply(2, function(x) ifelse(x > 0, 1, 0)) 
        my.dH <- ape::drop.tip(phy = host_tree, setdiff(host_tree$tip.label, phyloseq::sample_names(my.ps))) %>% ape::cophenetic.phylo()
        my.dP <- ape::cophenetic.phylo(phyloseq::phy_tree(my.ps))
        # run parafit
        my.pf <- try(ape::parafit(host.D = my.dH, para.D = my.dP, HP = my.HP, 
                                  nperm = nperm, correction = "cailliez", test.links = FALSE, silent = T), silent = T)
        if(class(my.pf) != "try-error"){
          # output results
          res.out[cluster, "stat"] <- my.pf$ParaFitGlobal 
          res.out[cluster, "P"] <- my.pf$p.global 
          #ps.out[[paste0("clust.",cluster)]] <- my.ps
        }  
      } 
    } 
  }
  res.out
}

##### Function to combine rare phyla and classes into "Other" for plotting
combn_rare_class <- function(df, min.p = 20, min.c = 10){
  # Options:
  # - df: data frame with taxonomy, and optionally other, columns including at least "Phylum" and "Class"
  # - min.p: minimum frequency of phylum
  # - min.c: minimum frequency of class
  #
  # initiate list to hold results for each phylum
  my.list <- list()
  # loop through Phyla
  for(P in sort(unique(df$Phylum))){
    # make phylum data frame
    my.df <- df[which(df$Phylum == P),]
    # replace any NA in class with "Unknown <phylum>"
    my.df[which(is.na(my.df$Class)), "Class"] <-  paste("Unknown", P)
    # if there are < min.p instances of Phylum, change to "Other"
    if(nrow(my.df) < min.p){
      my.df[,"Phylum"] <- "Other"
      my.df[,"Class"] <- "Other"
      my.list[[P]] <- my.df
    } else {
      # otherwise, loop through each class in the phylum
      for (C in sort(unique(my.df[,"Class"]))){
        # ignore unknowns
        if(C == paste("Unknown", P)){
          next
        } else {
          # if there are < min.c instances of the class, change to "Other <phylum>"
          if(length(which(my.df$Class == C)) < min.c){
            my.df[which(my.df$Class == C),"Class"] <- paste("Other", P)
          }
        }
      }
      # add to list
      my.list[[P]] <- my.df
    }
  }
  # combine list into one df
  df <- do.call(rbind, my.list)
  # sort df
  df <- df[with(df, order(Phylum, Class)), ]
  rownames(df) <- NULL
  # df
  df
}
