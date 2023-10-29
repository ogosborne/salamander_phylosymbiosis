library(phyloseq)
library(DESeq2)
library(dplyr)
load("results/4.Rarefy_and_subset/sal_unrar.RData")
# We cannot test the model: Species + Locality + Habitat + Bd because this model is not full rank.
# Instead, we use the model:  ~ Species + Site + Bd.
# make deseq2 object
dds <- phyloseq_to_deseq2(sal_f, ~ Species + Site + Bd)
# use poscounts to estimate size factors. Standard estimator fails due to zero-inflation.
dds = DESeq(dds, fitType="local", sfType = "poscounts")
# extract results
res = results(dds, contrast = c("Bd", "yes", "no"))
#dim(res2)
res = res[order(res$padj, na.last=T), ]
#dim(res)
res = cbind(as(res, "data.frame"), as(tax_table(sal_f)[rownames(res), ], "matrix"))
dim(res)
# add significance column
alpha = 0.05
res$signif <- res$padj <= alpha
res$signif_highnonBd <- res$padj <= alpha & res$log2FoldChange < 0
# load antiBd
antiBd <- read.csv("results/5.Antifungal_DB/blast.filt.uniq.csv", header = T)
# add antiBd column to 
res$antiBd <- res$ASV %in% antiBd$qseqid
# test association between known antiBd and differential abundance to Bd
mytab <- table(res[,c("antiBd", "signif")])
fish <- fisher.test(mytab)
# test with only those which are more abundant in uninfected
mytab2 <- table(res[,c("antiBd", "signif_highnonBd")])
fish2 <- fisher.test(mytab2)
# save results
dir.create("results/11.diff.abundance")
save(list = c("dds", "res", "fish", "fish2"), file = "results/11.diff.abundance/BdStatus_DESeq2.RData")
## make combined table of DESeq2 and Specificity results
load("results/11.diff.abundance/BdStatus_DESeq2.RData")
colnames(res)[c(1:6, 15:16)] <- paste0("DESeq2_",colnames(res)[c(1:6, 15:16)])
colnames(res)[17] <- "AntiBd"
spec_res <- read.csv("results/10.Specificity/Specificity_results.csv")
merged_res <- merge(res, spec_res, all = T)
# get seqs
seqs <- as(refseq(sal_f), "character")
seqsdf <- data.frame(ASV = names(seqs), Sequence = unname(seqs))
merged_res <- merge(merged_res, seqsdf, all = T)
# save table
write.csv(merged_res, "results/11.diff.abundance/DESeq_and_Spec_results.csv", row.names = F)


