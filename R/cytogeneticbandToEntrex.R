#' Map Cytogenetic Band to Entrez Gene IDs
#'
#' This function maps a specified cytogenetic band to corresponding Entrez gene IDs.
#'
#' @param cytogen_band The cytogenetic band to map.
#' @return A vector of Entrez gene IDs corresponding to the specified cytogenetic band.
#' @export
#'
#' @examples
#' # Map cytogenetic band to Entrez gene IDs
#' result <- cytogeneticbandToEntrex("Xq28")
#' head(result)
cytogeneticbandToEntrex <- function(cytogen_band) {

  # Gets the cytogene mapping based on the entrex for 10 genes.
  top_genes$chromosomal_band <- mapIds(org.Hs.eg.db,
                                       keys = top_genes$entrezid,
                                       column = "MAP",
                                       keytype = "ENTREZID",
                                       multiVals = "first")

  chromband_to_entrezid <- top_genes %>% filter(chromosomal_band == cytogen_band)
  entrezid_chromband <- chromband_to_entrezid$entrezid

  return(entrezid_chromband)
}
