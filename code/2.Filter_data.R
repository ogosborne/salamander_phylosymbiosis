library(dplyr)
library(phyloseq)
library(decontam)
library(Biostrings)
### CREATE PHYLOSEQ OBJECT ####
####### Import ASV table
featureTab <- read.csv("results/1.DADA2_and_taxonomy/feature_table.csv", header = T, row.names = 1)
featureTab <- otu_table(featureTab, taxa_are_rows = TRUE)
## Import Taxonomy
taxonomy <- as.matrix(read.csv("results/1.DADA2_and_taxonomy/taxonomy.csv", row.names = 1))
taxonomy <- tax_table(taxonomy)
## Import Metadata
metadata2020 <- read.csv("data/Salamanders_metadata.csv", header = T, row.names = 1)
# convert to phyloseq sample data
metadata2020 <- sample_data(metadata2020)
## Import Sequences
sequences <- readDNAStringSet("results/1.DADA2_and_taxonomy/asv_seqs.fasta")
# Make phyloseq object
all_f <- merge_phyloseq(featureTab, taxonomy, metadata2020, sequences) 
# nsamples(all_f) == 376, sum(sample_sums(all_f)) == 8489126
# Filter singletons 
all_f <- filter_taxa(all_f, function (x) {sum(x > 0) >1}, prune=TRUE) 
# nsamples(all_f) == 376, sum(sample_sums(all_f)) == 7792027 
## Remove NA's from Phylum level (might be artifacts) and Archeae
all_f <- subset_taxa(all_f, Phylum!="NA")
all_f <- subset_taxa(all_f, Phylum!="Cyanobacteria/Chloroplast")
all_f <- subset_taxa(all_f, Kingdom!="Archaea")
# nsamples(all_f) == 376, sum(sample_sums(all_f)) == 7287975
## Add ASVs as their own column in the tax_table
tax_table(all_f) <- cbind(tax_table(all_f), rownames(tax_table(all_f)))
## Rename taxa 
colnames(tax_table(all_f)) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "ASV")
# Filter out MOC community
all_f <- all_f %>% subset_samples( Site != "MOC")
all_f <- filter_taxa(all_f, function(x) sum(x) !=0, TRUE)
# nsamples(all_f) == 363, sum(sample_sums(all_f)) == 7250233
# Identify contaminants by prevalence
sample_data(all_f)$is.neg <- sample_data(all_f)$Sample_Type == "Control"
contamdf <- isContaminant(all_f, method="prevalence", neg="is.neg", threshold=0.1)
all_f <- prune_taxa(!contamdf$contaminant, all_f)
# nsamples(all_f) == 363, sum(sample_sums(all_f)) == 7206678
# Remove controls
all_f <- all_f %>%
  subset_samples(Sample_Type != "Control")
# nsamples(all_f) == 359, sum(sample_sums(all_f)) == 7205789
# Filter low abundant taxa 
all_f <- prune_taxa(taxa_sums(all_f) >= 10, all_f)
# nsamples(all_f) == 359, sum(sample_sums(all_f)) == 7193836
# write fasta of remaining ASVs to make tree
dir.create("results/3.ASV_tree/")
writeXStringSet(phyloseq::refseq(all_f), "results/3.ASV_tree/all_f.fasta")
# save phyloseq object
dir.create("results/2.Filter_data")
save(all_f, file = "results/2.Filter_data/filtered.RData")

