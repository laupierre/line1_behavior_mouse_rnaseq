library (gprofiler2)
library (openxlsx)

deg <- read.xlsx ("striatum_deseq2_tespex_OTvsFR_differential_expression.xlsx")
deg <- deg [deg$padj < 0.05, ]
deg_up <- deg[deg$log2FoldChange > 0, ]
deg_up <- deg_up$gene_name
head (deg_up)


gp_up <- gost(query = deg_up, organism = "mmusculus", numeric_ns = "ENTREZGENE_ACC", sources = c("GO:BP"), exclude_iea=TRUE, evcodes = TRUE, significant = FALSE)
gp_up <- gp_up$result[ ,c("query", "significant", "p_value", "term_name", "term_size", "intersection_size", "intersection")]
gp_up$query <- "up"
head(gp_up)

