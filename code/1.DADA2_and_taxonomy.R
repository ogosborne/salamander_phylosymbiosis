# load libraries
library(dada2)
library(ggplot2)
library(gridExtra)
library(phyloseq)
library(Biostrings)
# Directories with the 2 different runs
rundirs <- c(DS1 = "data/2020_data/DS1", DS2 = "data/2020_data/DS2")
runs <- names(rundirs)
# Create output directory
dir.create("results/1.DADA2_and_taxonomy/")
# Get file names
files <- list()
for(i in runs){
  files[[i]] <- list()
  for(r in 1:2){
    file.names <- sort(list.files(rundirs[i], pattern = paste0("_R", r, "_001.fastq.gz"), full.names = TRUE))
    samp.names <- sapply(strsplit(basename(file.names), "_"), `[`, 1)
    filt.names <- file.path(rundirs[i], "filtered", paste0(samp.names, "_filt_R", r,".fastq.gz"))
    my.list <- list(file.names = file.names, samp.names = samp.names, filt.names = filt.names)
    files[[i]][[paste0("R",r)]] <- my.list
  }
}
# Plot quality profiles to decide trimming cut-offs
p1 <- plotQualityProfile(files$DS1$R1$file.names, aggregate = T)
p2 <- plotQualityProfile(files$DS1$R2$file.names, aggregate = T)
p3 <- plotQualityProfile(files$DS2$R1$file.names, aggregate = T)
p4 <- plotQualityProfile(files$DS2$R2$file.names, aggregate = T)
p <- grid.arrange(p1 + ggtitle("DS1 R1"), 
                        p2 + ggtitle("DS1 R2"),
                        p3 + ggtitle("DS2 R1"),
                        p4 + ggtitle("DS2 R2"))
ggsave(plot = p, file = "results/1.DADA2_and_taxonomy/quality_profiles.pdf", device = "pdf")
# Filter and trim
trunclen <- list(DS1 = c(270,180), DS2 = c(270,180)) 
trimlog <- list()
for(i in runs){
  trimlog[[i]] <-  filterAndTrim(files[[i]]$R1$file.names, files[[i]]$R1$filt.names, files[[i]]$R2$file.names, files[[i]]$R2$filt.names, 
                                 truncLen=trunclen[[i]],
                                 maxN=0, 
                                 maxEE=c(2,2), 
                                 trimLeft = 19, 
                                 trimRight = 23,
                                 truncQ=2, 
                                 rm.phix=TRUE,
                                 compress=TRUE, 
                                 multithread=TRUE)
}
save(trimlog, file = "results/1.DADA2_and_taxonomy/trimlog.RData")
# write read counts
log.comb <- data.frame(trimlog$DS1 + trimlog$DS2)
log.comb$samp <- row.names(log.comb)
log.comb <- log.comb[,c("samp","reads.in","reads.out")]
write.csv(DS.comb, file = "results/1.DADA2_and_taxonomy/readcounts.csv", row.names = F)
# Learn the Error Rates
err <- list()
for(i in runs){
  err[[i]] <- list()
  for(r in 1:2){
    err[[i]][[paste0("R", r)]] <- learnErrors(files[[i]][[paste0("R",r)]]$filt.names, multithread=TRUE)
  }
}
save(err, file = "results/1.DADA2_and_taxonomy/err.RData")
# derep
derep <- list()
for(i in runs){
  derep[[i]] <- list()
  for(r in 1:2){
    my.derep <- derepFastq(files[[i]][[paste0("R",r)]]$filt.names)
    names(my.derep) <- files[[i]][[paste0("R",r)]]$samp.names
    derep[[i]][[paste0("R", r)]] <- my.derep
  }
}
save(derep, file = "results/1.DADA2_and_taxonomy/derep.RData")
# Sample Inference
dadares <- list()
for(i in runs){
  dadares[[i]] <- list()
  for(r in 1:2){
    my.dadares <- dada(derep[[i]][[paste0("R",r)]], err = err[[i]][[paste0("R", r)]], multithread = TRUE)
    dadares[[i]][[paste0("R", r)]] <- my.dadares
  }
}
save(dadares, file = "results/1.DADA2_and_taxonomy/dadares.RData")
# Merge paired reads
merged <- list()
for(i in runs){
  dada1 <- dadares[[i]][["R1"]]
  dada2 <- dadares[[i]][["R2"]]
  derep1 <- derep[[i]][["R1"]]
  derep2 <- derep[[i]][["R2"]]
  my.merged <- mergePairs(dada1, derep1, dada2, derep2, verbose=TRUE)
  merged[[i]] <- my.merged
}
save(merged, file = "results/1.DADA2_and_taxonomy/merged.RData")
# Construct sequence table
seqtab <- list()
for(i in runs){
  seqtab[[i]] <- makeSequenceTable(merged[[i]])
}
save(seqtab, file = "results/1.DADA2_and_taxonomy/seqtab.RData")
# Merge multiple runs
seqtab.all <- mergeSequenceTables(seqtab$DS1, seqtab$DS2, repeats = "sum")
save(seqtab.all, file = "results/1.DADA2_and_taxonomy/seqtab.all.RData")
# Remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab.all, method="consensus", multithread=TRUE, verbose=TRUE)
save(seqtab.nochim, file = "results/1.DADA2_and_taxonomy/seqtab.nochim.RData")
# Assign Taxonomy
taxa <- assignTaxonomy(seqs = seqtab.nochim, refFasta = "data/RDP_train_set/rdp_train_set_16.fa.gz", multithread=TRUE, verbose = T)
save(taxa, file = "results/1.DADA2_and_taxonomy/taxa1.RData")
# Assign species
chunk.size <- 5000  # size of taxonomy increments
#https://github.com/chuvpne/dada2-pipeline/commit/7964cd67da52faadd31f8d93da6385741984360b
taxonomy.species <- do.call(rbind,
                            lapply(split(c(1:nrow(taxa)), sort(c(1:nrow(taxa))%%ceiling(nrow(taxa)/chunk.size))), function(x){return(addSpecies(taxa[x, ], "data/RDP_train_set/rdp_species_assignment_16.fa.gz"))}))
save(taxonomy.species, file = "results/1.DADA2_and_taxonomy/taxa2.RData")
#  make phyloseq object to combine feature table and taxonomy table in same order
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
         tax_table(taxonomy.species))
# rename samples to remove '-'
sample_names(ps) <- gsub("-","_",sample_names(ps))
## extract seqs and rename ASVs 
new.names <- paste0("ASV", seq(ntaxa(ps))) # Define new names ASV1, ASV2, ...
seqs <- taxa_names(ps) # Store sequences
names(seqs) <- new.names # Name sequences
seqs <- DNAStringSet(seqs, use.names = T) # convert to DNAStringSet
taxa_names(ps) <- new.names # Rename ASVs
## extract feature table
ftab <- t(as(otu_table(ps), "matrix"))
## extract tax table
tax <- as(tax_table(ps), "matrix")
## write to files
write.csv(ftab, "results/1.DADA2_and_taxonomy/feature_table.csv")
write.csv(tax, "results/1.DADA2_and_taxonomy/taxonomy.csv")
writeXStringSet(seqs, "results/1.DADA2_and_taxonomy/asv_seqs.fasta")
