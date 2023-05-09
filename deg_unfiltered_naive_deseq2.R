## 210825_A00558_0141_AH232YDMXY

library (openxlsx)
library (DESeq2)
library (ggplot2)
library (ggrepel)
library (pheatmap)


anno <- read.delim ("gencode.vM32.annotation.txt")
anno <- anno[ ,grep ("transcript_id", colnames (anno), invert=TRUE)]
anno <- unique (anno)

a <- read.delim ("subread.counts.txt", skip=1)
a <- a[ ,grep ("ene|bam", colnames (a))]
a <- a[grep ("miRNA|Mt_tRNA|Mt_rRNA|rRNA|snRNA|snoRNA|scRNA|sRNA|misc_RNA|scaRNA|ribozyme|IG_|TR_", a$gene_type, invert=TRUE), ]
colnames (a) <- gsub ("_S[0-9]+.*", "", colnames (a)) 
colnames (a) <- gsub (".*TRM_", "", colnames (a))

a <- merge (a, anno, by.x="Geneid", by.y="gene_id", all.x=TRUE) 
a <- a[ ,grep ("gene_type.y|gene_name.y", colnames (a), invert=TRUE)]
colnames (a) [colnames (a) == "gene_name.x"] <- "gene_name"
colnames (a) [colnames (a) == "gene_type.x"] <- "gene_type"

#write.xlsx (a, "star_gene_raw_counts.xlsx", rowNames=F)


annot <- a
annot <- annot[ ,c("Geneid", "gene_name", "gene_type", "mgi_id", "external_gene_name", "description")]


torm <- c("gene_name", "gene_type", "mgi_id", "external_gene_name", "description")
a <- a[ ,!colnames (a) %in% torm]
row.names (a) <- a[ ,1]
colnames (a) <- gsub ("star.IIT_", "", colnames (a))
a <- a[ ,-1]



## Naive experiment

#FR: 1623, 1636, 1637, 1638, 1971, 1972, 1973
#ST: 1627, 1628, 1629, 1631, 1962, 1966, 1950
#OR: 1616, 1617, 1618, 1622, 1947, 1949, 1952

pheno <- data.frame (sample=c("RNT_1623", "RNT_1636", "RNT_1637", "RNT_1638", "RNT_1971", "RNT_1972", "RNT_1973",
            "RNT_1627", "RNT_1628", "RNT_1629", "RNT_1631", "RNT_1962", "RNT_1966", "RNT_1950",
            "RNT_1616", "RNT_1617", "RNT_1618", "RNT_1622", "RNT_1947", "RNT_1949", "RNT_1952"))
            
pheno$genotype <- c(rep ("FR", 7), rep ("ST", 7), rep ("OR", 7))


a <- a[ ,colnames (a) %in% pheno$sample]

idx <- match (colnames (a), pheno$sample)
pheno <- pheno[idx, ]

stopifnot (colnames (a) == pheno$sample)



counts <- a

dds <- DESeqDataSetFromMatrix(countData = round (counts), colData = pheno, design = ~ genotype)
                                 
# keep <- rowSums(counts(dds)) >= 50
keep <- rowSums(counts(dds) >= 50) >= dim (counts)[2]/2
dds <- dds[keep,]
dds

# first contrast of interest (ST vs FR)
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, contrast=list("genotype_ST_vs_FR"))

res <- merge (data.frame (res), counts (dds), by="row.names")
#res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

write.xlsx (res, "striatum_deseq2_STvsFR_differential_expression.xlsx", rowNames=F)

boxplot (res$log2FoldChange)
abline (h=0)



# second contrast of interest (OR vs FR)

res <- results(dds, contrast=list("genotype_OR_vs_FR"))

res <- merge (data.frame (res), counts (dds), by="row.names")
#res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid")
colnames (res)[1] <- "Geneid"
res <- res[order (res$padj), ]

write.xlsx (res, "striatum_deseq2_ORvsFR_differential_expression.xlsx", rowNames=F)

boxplot (res$log2FoldChange)
abline (h=0)










