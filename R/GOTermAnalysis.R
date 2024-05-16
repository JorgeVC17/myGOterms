#' GO Term Enrichment Analysis
#'
#' Perform GO term enrichment analysis based on differential gene expression results.
#'
#' @param dge_results A data frame containing differential gene expression results.
#' @param upregulated Logical indicating whether to analyze upregulated genes.
#' @param lfc_value Log fold change threshold for differential expression.
#' @param padj_value Adjusted p-value threshold for differential expression.
#' @return A summary of the GO term enrichment analysis.
#' @export
#'
#' @examples
#' # Perform GO term enrichment analysis on upregulated genes
#' result <- GOTermAnalysis(dge_results, upregulated = TRUE, lfc_value = 1, padj_value = 0.05)

GOTermAnalysis <- function(dge_results, upregulated, lfc_value, padj_value) {
  # Make a list of all genes in the dataset
  all_genes <- dge_results %>% data.frame() %>% rownames()
  if (upregulated==TRUE) {
    # Make a list of all upregulated genes in the dataset
    upregulated_genes <- dge_results %>% data.frame() %>%
      filter(log2FoldChange > lfc_value, padj < padj_value) %>% rownames()
    # Execute GO term enrichment analysis
    GO_terms_up <- new("GOHyperGParams",
                                       geneIds = upregulated_genes,
                                       universeGeneIds = all_genes,
                                       annotation = "org.Hs.eg.db",
                                       ontology = "BP",
                                       pvalueCutoff = 1,
                                       testDirection = "over")

    goterm_analysis_up_summary <- summary(hyperGTest(GO_terms_up))
    return(goterm_analysis_up_summary)
  } else {
    # Make a list of all downregulated genes in the dataset
    downregulated_genes <- dge_results %>% data.frame() %>%
      filter(log2FoldChange < -lfc_value, padj < padj_value) %>% rownames()
    #Uitvoer GO term enrichment analysis
    GO_terms_down <- new("GOHyperGParams",
                                         geneIds = downregulated_genes,
                                         universeGeneIds = all_genes,
                                         annotation = "org.Hs.eg.db",
                                         ontology = "BP",
                                         pvalueCutoff = 1,
                                         testDirection = "over")

    goterm_analysis_down_summary <- summary(hyperGTest(GO_terms_down))
    return(goterm_analysis_down_summary)
  }
}
