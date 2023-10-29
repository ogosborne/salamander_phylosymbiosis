library(ape)
library(phyloseq)
library(Biostrings)

load("results/2.Filter_data/filtered.RData")
# import trees and add to phyloseq objects, multi2di resolves any polytomies
phy_tree(all_f) <- multi2di(read_tree("results/3.ASV_tree/all_f.nwk"))
# nsamples(all_f) == 359 , sum(sample_sums(all_f)) == 7193836
# remove low count samples
all_f <- prune_samples(sample_sums(all_f ) > 2900, all_f )
# nsamples(all_f) == 340, sum(sample_sums(all_f)) == 7167421
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
