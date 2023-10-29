library(rBLAST)
library(Biostrings)
# load sequences
seq <- readDNAStringSet("results/5.Antifungal_DB/all_f.fasta")
# make blast DB
makeblastdb(file = "data/AmphibBacJune15.2020_All_Inhibitory.fasta", dbtype = "nucl", args = "-parse_seqids -out results/5.Antifungal_DB/woodhams")
# link to blast DB
bl <- blast(db="results/5.Antifungal_DB/woodhams")
# run blast
result <- predict(bl, seq, 
              BLAST_args = "-perc_identity 100 -num_threads 6",
              custom_format = "qseqid sseqid qlen slen qstart qend sstart send length pident qcovs evalue bitscore") 
# filter by query coverage
result.filt <- result[which(result$qcovs == 100),]
# get one hit per query
result.filt.uniq <- result.filt[!duplicated(result.filt$qseqid),]
# output results
write.csv(result, file = "results/5.Antifungal_DB/blast.raw.csv", row.names = F)
write.csv(result.filt, file = "results/5.Antifungal_DB/blast.filt.csv", row.names = F)
write.csv(result.filt.uniq, file = "results/5.Antifungal_DB/blast.filt.uniq.csv", row.names = F)