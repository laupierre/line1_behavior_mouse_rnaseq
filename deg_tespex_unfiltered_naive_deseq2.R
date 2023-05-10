## 210825_A00558_0141_AH232YDMXY

library (openxlsx)
library (DESeq2)
library (ggplot2)
library (ggrepel)
library (pheatmap)
library (dplyr)

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


## TEspeX result

tesp <- read.delim ("outfile.txt", row.names=1)
tesp$transname <- gsub ("#.*", "", row.names (tesp))
tesp$transname <- transname <- gsub ("_[3|5|o].*", "", tesp$transname)

transfamily <- gsub (".*/", "", row.names (tesp))
annot_trans <- unique (data.frame (cbind (transname, transfamily)))


# aggregate multiple columns
tesp <- data.frame (tesp %>% group_by (transname) %>% summarise(across(everything(), sum)))

# same as in base: tesp <- aggregate(. ~ transname, tesp, sum)
# data.frame (t (tesp[tesp$transname == "L1MdA_II", ]))

colnames (tesp) <- gsub ("IIT_", "", colnames (tesp))
colnames (tesp) <- gsub ("_S.*", "", colnames (tesp))
row.names (tesp) <- tesp[ ,1]
tesp <- tesp[ ,-1]

tesp <- tesp[ ,colnames (tesp) %in% pheno$sample]

idx <- match (colnames (tesp), pheno$sample)
tesp <- tesp[ , idx]

stopifnot (colnames (tesp) == pheno$sample)




counts <- rbind (a, tesp)

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
res <- merge (res, annot, by.x="Row.names", by.y="Geneid", all.x=TRUE)
colnames (res)[1] <- "Geneid"
res <- resa <- res[order (res$padj), ]

idx2 <- which (is.na (res$external_gene_name))
res$external_gene_name [idx2] <- res$Geneid[is.na (res$external_gene_name)]
res$gene_name [idx2] <- res$Geneid[idx2]
res$gene_type [idx2] <- "transposon"

for (i in (1:dim (res)[1])) {
if (res$gene_type[i] == "transposon") {
print (res$gene_name[i])
#res$description[i] <- annot_trans$transfamily [annot_trans$transname == res$gene_name[i]]
res$description[i]  <- annot_trans$transfamily [annot_trans$transname == res$gene_name[i]]
}
}

write.xlsx (res, "striatum_deseq2_tespex_STvsFR_differential_expression.xlsx", rowNames=F)

boxplot (res$log2FoldChange)
abline (h=0)




# second contrast of interest (OR vs FR)

res <- results(dds, contrast=list("genotype_OR_vs_FR"))

res <- merge (data.frame (res), counts (dds), by="row.names")
#res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid", all.x=TRUE)
colnames (res)[1] <- "Geneid"
res <- resa <- res[order (res$padj), ]

idx2 <- which (is.na (res$external_gene_name))
res$external_gene_name [idx2] <- res$Geneid[is.na (res$external_gene_name)]
res$gene_name [idx2] <- res$Geneid[idx2]
res$gene_type [idx2] <- "transposon"

for (i in (1:dim (res)[1])) {
if (res$gene_type[i] == "transposon") {
print (res$gene_name[i])
#res$description[i] <- annot_trans$transfamily [annot_trans$transname == res$gene_name[i]]
res$description[i]  <- annot_trans$transfamily [annot_trans$transname == res$gene_name[i]]
}
}

write.xlsx (res, "striatum_deseq2_tespex_ORvsFR_differential_expression.xlsx", rowNames=F)

boxplot (res$log2FoldChange)
abline (h=0)




## PCA plot

vsd <- vst(dds, blind=FALSE)
pcaData <- plotPCA(vsd, intgroup=c("genotype"), returnData=TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))

ggplot(pcaData, aes(PC1, PC2, color=genotype, shape=genotype)) +
  		geom_point(size=3) +
  		xlab(paste0("PC1: ",percentVar[1],"% variance")) +
  		ylab(paste0("PC2: ",percentVar[2],"% variance")) + 
		coord_fixed ()

ggsave ("PCA plot naive experiment.pdf")



## Comparison with previous WGCNA based pipeline

prev <- read.xlsx ("/Volumes/texas/iit_projects/tonini_version2/striatum_deseq2_STvsFR_differential_expression_wgcna_pipeline.xlsx")
prev <- merge (resa, prev, by.x="gene_name", by.y="Geneid")
plot (prev$log2FoldChange.x, prev$log2FoldChange.y, col=ifelse (prev$padj.x < 0.05 & prev$padj.y < 0.05, "darkblue", "black"), xlab="new_pipe_limma", ylab="prev_pipe_limma")
abline (h=0)
abline (v=0)
abline (0,1)

cor (prev$log2FoldChange.x, prev$log2FoldChange.y, method="pearson")
# 0.99


prev <- read.xlsx ("/Volumes/texas/iit_projects/tonini_version2/striatum_deseq2_ORvsFR_differential_expression_wgcna_pipeline.xlsx")
prev <- merge (resb, prev, by.x="gene_name", by.y="Geneid")
plot (prev$log2FoldChange.x, prev$log2FoldChange.y, col=ifelse (prev$padj.x < 0.05 & prev$padj.y < 0.05, "darkblue", "black"), xlab="new_pipe_limma", ylab="prev_pipe_limma")
abline (h=0)
abline (v=0)
abline (0,1)

cor (prev$log2FoldChange.x, prev$log2FoldChange.y, method="pearson")
# 0.99


