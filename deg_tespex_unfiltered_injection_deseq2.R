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



## Injection experiment

#ST with shCTRL: 392, 393, 394, 395, 396, 397, 398
#OT with shCTRL: 362, 363, 364, 365, 368, 370, 371
#OT with shL1a:  342, 343, 344, 346, 359, 360, 361

pheno <- data.frame (sample=c("RNT_392", "RNT_393", "RNT_394", "RNT_395", "RNT_396", "RNT_397", "RNT_398",
            "RNT_362", "RNT_363", "RNT_364", "RNT_365", "RNT_368", "RNT_370", "RNT_371",
            "RNT_342", "RNT_343", "RNT_344", "RNT_346", "RNT_359", "RNT_360", "RNT_361"))
            
pheno$genotype <- c(rep ("STCTRL", 7), rep ("OTCTRL", 7), rep ("OTL1", 7))


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


# first contrast of interest (OTL1 vs OTCTRL)
dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds, contrast=list("genotype_OTL1_vs_OTCTRL"))

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
#print (res$gene_name[i])
res$description[i]  <- annot_trans$transfamily [annot_trans$transname == res$gene_name[i]]
}
}

write.xlsx (res, "striatum_deseq2_tespex_OTSHL1vsOTSHCTRL_injection_differential_expression.xlsx", rowNames=F)

boxplot (res$log2FoldChange)
abline (h=0)



# second contrast of interest (STCTRL_vs_OTCTRL)

rm (res)
res <- results(dds, contrast=list("genotype_STCTRL_vs_OTCTRL"))

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
#print (res$gene_name[i])
res$description[i]  <- annot_trans$transfamily [annot_trans$transname == res$gene_name[i]]
}
}

write.xlsx (res, "striatum_deseq2_tespex_STSHCTRLvsOTSHCTRL_injection_differential_expression.xlsx", rowNames=F)

boxplot (res$log2FoldChange)
abline (h=0)



# third contrast from the two previous contrats, using the list() style of contrasts !!!

res <- results(dds, contrast=list("genotype_OTL1_vs_OTCTRL", "genotype_STCTRL_vs_OTCTRL"))     ## this is equivalent to OTL1/STCTRL

res <- merge (data.frame (res), counts (dds), by="row.names")
#res <- merge (data.frame (res), round (counts (dds, normalized=TRUE)), by="row.names")
res <- merge (res, annot, by.x="Row.names", by.y="Geneid", all.x=TRUE)
colnames (res)[1] <- "Geneid"
res <- resb <- res[order (res$padj), ]

idx2 <- which (is.na (res$external_gene_name))
res$external_gene_name [idx2] <- res$Geneid[is.na (res$external_gene_name)]
res$gene_name [idx2] <- res$Geneid[idx2]
res$gene_type [idx2] <- "transposon"

for (i in (1:dim (res)[1])) {
if (res$gene_type[i] == "transposon") {
#print (res$gene_name[i])
res$description[i]  <- annot_trans$transfamily [annot_trans$transname == res$gene_name[i]]
}
}

write.xlsx (res, "striatum_deseq2_tespex_OTSHL1vsSTSHCTRL_injection_differential_expression.xlsx", rowNames=F)

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

ggsave ("PCA plot injection experiment.pdf")











