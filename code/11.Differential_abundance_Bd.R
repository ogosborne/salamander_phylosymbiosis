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
#dim(res2)
res = res[order(res$padj, na.last=T), ]
#dim(res)
res = cbind(as(res, "data.frame"), as(tax_table(sal_f)[rownames(res), ], "matrix"))
dim(res)
# add significance column
alpha = 0.05
res$signif <- res$padj <= alpha
res[which(is.na(res$signif)),"signif"] <- FALSE
res$signif_highnonBd <- res$signif == TRUE & res$log2FoldChange < 0
res$signif_lownonBd <- res$signif == TRUE & res$log2FoldChange > 0
# load antiBd
antiBd <- read.csv("results/5.Antifungal_DB/blast.filt.uniq.csv", header = T)
# add antiBd column to 
res$antiBd <- res$ASV %in% antiBd$qseqid
# test association between known antiBd and differential abundance to Bd
mytab <- table(res[,c("signif", "antiBd")])
apply(mytab, MARGIN = 2, function(x) x[2]/sum(x) * 100)
fish <- fisher.test(mytab) ; fish
# test with only those which are more abundant in uninfected
mytab2 <- table(res[,c("signif_highnonBd", "antiBd")])
apply(mytab2, MARGIN = 2, function(x) x[2]/sum(x) * 100)
fish2 <- fisher.test(mytab2) ; fish2
# test with only those which are more abundant in infected
mytab3 <- table(res[,c("signif_lownonBd", "antiBd")])
apply(mytab3, MARGIN = 2, function(x) x[2]/sum(x) * 100)
fish3 <- fisher.test(mytab3) ; fish3
# save results
dir.create("results/11.diff.abundance")
save(list = c("dds", "res", "fish", "fish2", "fish3"), file = "results/11.diff.abundance/BdStatus_DESeq2.RData")
## make combined table of DESeq2 and Specificity results
load("results/11.diff.abundance/BdStatus_DESeq2.RData")
colnames(res)[c(1:6, 15:17)] <- paste0("DESeq2_",colnames(res)[c(1:6, 15:17)])
colnames(res)[18] <- "AntiBd"
spec_res <- read.csv("results/10.Specificity/Specificity_results.csv")
merged_res <- merge(res, spec_res, all = T)
# get seqs
seqs <- as(refseq(sal_f), "character")
seqsdf <- data.frame(ASV = names(seqs), Sequence = unname(seqs))
merged_res <- merge(merged_res, seqsdf, all = T)
# save table
write.csv(merged_res, "results/11.diff.abundance/DESeq_and_Spec_results.csv", row.names = F)


