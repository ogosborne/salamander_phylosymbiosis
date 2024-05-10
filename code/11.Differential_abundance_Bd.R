library(phyloseq)
library(DESeq2)
library(dplyr)
load("results/4.Rarefy_and_subset/sal_unrar.RData")
# We cannot test the model: Species + Habitat + Locality + Bd because this model is not full rank.
# Instead, we combine Species and Habitat into one variable and use the formula ~ Species_Habitat + Locality + Bd.
# prepare explanatory factors
sample_data(sal_f)$Species_Habitat <- as.factor(gsub(" ", "", paste(sample_data(sal_f)$Species, sample_data(sal_f)$Habitat, sep = "_")))
sample_data(sal_f)$Locality <- as.factor(gsub(" ", "", sample_data(sal_f)$Locality))
sample_data(sal_f)$Bd <- as.factor(sample_data(sal_f)$Bd)
# make deseq2 object
dds <- phyloseq_to_deseq2(sal_f, ~ Species_Habitat + Locality + Bd)
# use poscounts to estimate size factors. Standard estimator fails due to zero-inflation.
dds = DESeq(dds, fitType="local", sfType = "poscounts")
# extract results
res = results(dds, contrast = c("Bd", "yes", "no"))
res = res[order(res$padj, na.last=T), ]
res = cbind(as(res, "data.frame"), as(tax_table(sal_f)[rownames(res), ], "matrix"))
# add significance column
alpha = 0.05
res$signif <- res$padj <= alpha
res[which(is.na(res$signif)),"signif"] <- FALSE
res$signif_highnonBd <- res$signif == TRUE & res$log2FoldChange < 0
res$signif_lownonBd <- res$signif == TRUE & res$log2FoldChange > 0
# load antiBd
antiBd <- read.csv("results/5.Antifungal_DB/blast.filt.uniq.csv", header = T)
# add antiBd column to 
res$AntiBd <- res$ASV %in% antiBd$qseqid
## make combined table of DESeq2 and Specificity results
colnames(res)[c(1:6, 15:17)] <- paste0("DESeq2_",colnames(res)[c(1:6, 15:17)])
spec_res <- read.csv("results/10.Specificity/Specificity_results.csv")
merged_res <- merge(res, spec_res, all = T)
# get seqs
seqs <- as(refseq(sal_f), "character")
seqsdf <- data.frame(ASV = names(seqs), Sequence = unname(seqs))
merged_res <- merge(merged_res, seqsdf, all = T)
# save table
write.csv(merged_res, "results/11.diff.abundance/DESeq_and_Spec_results.csv", row.names = F)
# test association between known antiBd and differential abundance to Bd
tab.abd.dabd <- table(merged_res[,c("DESeq2_signif", "AntiBd")])
apply(tab.abd.dabd, MARGIN = 2, function(x) x[2]/sum(x) * 100)
fish.abd.dabd <- fisher.test(tab.abd.dabd) ; fish.abd.dabd
# test with only those which are more abundant in uninfected
tab.abd.lobd <- table(merged_res[,c("DESeq2_signif_highnonBd", "AntiBd")])
apply(tab.abd.lobd, MARGIN = 2, function(x) x[2]/sum(x) * 100)
fish.abd.lobd <- fisher.test(tab.abd.lobd) ; fish.abd.lobd
# test with only those which are more abundant in infected
tab.abd.hibd <- table(merged_res[,c("DESeq2_signif_lownonBd", "AntiBd")])
apply(tab.abd.hibd, MARGIN = 2, function(x) x[2]/sum(x) * 100)
fish.abd.hibd <- fisher.test(tab.abd.hibd) ; fish.abd.hibd
# test association between known antiBd and significant specificity to Bd load (from specificity analysis)
BdSpecDF <- data.frame(sigSpecBdLoad =  merged_res$bdlo_dist.Pval < 0.05 & merged_res$bdlo_dist.Spec < 0, AntiBd = merged_res$AntiBd)
BdSpecDF <- BdSpecDF[complete.cases(BdSpecDF),]
tab.abd.spbd <- table(BdSpecDF)
apply(tab.abd.spbd, MARGIN = 2, function(x) x[2]/sum(x) * 100)
fish.abd.spbd <- fisher.test(tab.abd.spbd) ; fish.abd.spbd
# save results
dir.create("results/11.diff.abundance")
save(list = c("dds", "merged_res", "fish.abd.dabd", "fish.abd.lobd", "fish.abd.hibd", "fish.abd.spbd"), file = "results/11.diff.abundance/BdStatus_DESeq2.RData")

