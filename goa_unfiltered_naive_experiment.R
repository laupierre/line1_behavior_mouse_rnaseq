## See https://cran.r-project.org/web/packages/gprofiler2/vignettes/gprofiler2.html

library (gprofiler2)
library (openxlsx)


## first contrast

deg <- read.xlsx ("striatum_deseq2_tespex_OTvsFR_differential_expression.xlsx")
deg <- deg[deg$gene_type != "transposon", ]

custom_bg <- deg$gene_name

deg <- deg [deg$padj < 0.05, ]
deg_up <- deg[deg$log2FoldChange > 0, ]
deg_up <- deg_up$gene_name
head (deg_up)

gp_up <- gost(query = deg_up, organism = "mmusculus", numeric_ns = "ENTREZGENE_ACC", sources = c("GO:BP"), exclude_iea=TRUE, evcodes = TRUE, significant = TRUE,
              custom_bg = custom_bg, domain_scope = "custom_annotated")
gp_up <- gp_up$result[ ,c("term_id", "term_name", "significant", "p_value", "term_size", "query_size", "intersection_size", "parents", "intersection")]
write.xlsx (gp_up, "striatum_deseq2_tespex_OTvsFR_gene_ontology_analysis_up_genes.xlsx")


deg_down <- deg[deg$log2FoldChange < 0, ]
deg_down <- deg_down$gene_name
head (deg_down)

gp_down <- gost(query = deg_down, organism = "mmusculus", numeric_ns = "ENTREZGENE_ACC", sources = c("GO:BP"), exclude_iea=TRUE, evcodes = TRUE, significant = TRUE,
              custom_bg = custom_bg, domain_scope = "custom_annotated")
gp_down <- gp_down$result[ ,c("term_id", "term_name", "significant", "p_value", "term_size", "query_size", "intersection_size", "parents", "intersection")]
write.xlsx (gp_down, "striatum_deseq2_tespex_OTvsFR_gene_ontology_analysis_down_genes.xlsx")




## second contrast

rm (deg, gp_up, gp_down)

deg <- read.xlsx ("striatum_deseq2_tespex_STvsFR_differential_expression.xlsx")
deg <- deg[deg$gene_type != "transposon", ]

custom_bg <- deg$gene_name

deg <- deg [deg$padj < 0.05, ]
deg_up <- deg[deg$log2FoldChange > 0, ]
deg_up <- deg_up$gene_name
head (deg_up)

gp_up <- gost(query = deg_up, organism = "mmusculus", numeric_ns = "ENTREZGENE_ACC", sources = c("GO:BP"), exclude_iea=TRUE, evcodes = TRUE, significant = TRUE,
              custom_bg = custom_bg, domain_scope = "custom_annotated")
gp_up <- gp_up$result[ ,c("term_id", "term_name", "significant", "p_value", "term_size", "query_size", "intersection_size", "parents", "intersection")]
write.xlsx (gp_up, "striatum_deseq2_tespex_STvsFR_gene_ontology_analysis_up_genes.xlsx")


deg_down <- deg[deg$log2FoldChange < 0, ]
deg_down <- deg_down$gene_name
head (deg_down)

gp_down <- gost(query = deg_down, organism = "mmusculus", numeric_ns = "ENTREZGENE_ACC", sources = c("GO:BP"), exclude_iea=TRUE, evcodes = TRUE, significant = TRUE,
              custom_bg = custom_bg, domain_scope = "custom_annotated")
gp_down <- gp_down$result[ ,c("term_id", "term_name", "significant", "p_value", "term_size", "query_size", "intersection_size", "parents", "intersection")]
write.xlsx (gp_down, "striatum_deseq2_tespex_STvsFR_gene_ontology_analysis_down_genes.xlsx")




## third contrast

rm (deg, gp_up, gp_down)


deg <- read.xlsx ("striatum_deseq2_tespex_OTvsST_differential_expression.xlsx")
deg <- deg[deg$gene_type != "transposon", ]

custom_bg <- deg$gene_name

deg <- deg [deg$padj < 0.05, ]
deg_up <- deg[deg$log2FoldChange > 0, ]
deg_up <- deg_up$gene_name
head (deg_up)

gp_up <- gost(query = deg_up, organism = "mmusculus", numeric_ns = "ENTREZGENE_ACC", sources = c("GO:BP"), exclude_iea=TRUE, evcodes = TRUE, significant = TRUE,
              custom_bg = custom_bg, domain_scope = "custom_annotated")
gp_up <- gp_up$result[ ,c("term_id", "term_name", "significant", "p_value", "term_size", "query_size", "intersection_size", "parents", "intersection")]
write.xlsx (gp_up, "striatum_deseq2_tespex_OTvsST_gene_ontology_analysis_up_genes.xlsx")


deg_down <- deg[deg$log2FoldChange < 0, ]
deg_down <- deg_down$gene_name
head (deg_down)

gp_down <- gost(query = deg_down, organism = "mmusculus", numeric_ns = "ENTREZGENE_ACC", sources = c("GO:BP"), exclude_iea=TRUE, evcodes = TRUE, significant = TRUE,
              custom_bg = custom_bg, domain_scope = "custom_annotated")
gp_down <- gp_down$result[ ,c("term_id", "term_name", "significant", "p_value", "term_size", "query_size", "intersection_size", "parents", "intersection")]
write.xlsx (gp_down, "striatum_deseq2_tespex_OTvsST_gene_ontology_analysis_down_genes.xlsx")








