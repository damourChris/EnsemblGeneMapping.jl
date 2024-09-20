#' Extracts the first gene symbol from a given column
#' in an ExpressionSet object.
#'
#' This function takes an ExpressionSet object (`eset`) and a column name
#' (`gene_col`), and extracts the first gene symbol from each entry in
#' the specified column. The gene symbols are separated by a specified
#' separator (`separator`). The function returns a new ExpressionSet
#' with the original column replaced with the  extracted  gene symbols.
#'
#' @param eset An ExpressionSet object.
#' @param gene_col The name of the column containing gene symbols.
#' @param separator The separator used to separate multiple gene symbols
#' in each entry. Default is '///'
#'
#' @return A new ExpressionSet object with the first gene symbol
#' extracted from each entry.
#'
#' @examples
#' # Extract the first gene symbol from the "genes" column
#' eset_modified <- extract_first_gene_symbol(eset, "Entrez_gene_id", "///")
#'
#' # Print the modified ExpressionSet object
#' print(eset_modified)
#'
#' @importfrom Biobase fData
extract_first_gene_symbol <- function(eset, gene_col, separator = "///") {
  base_eset <- eset

  genes_to_map_raw <- fData(base_eset)[[gene_col]]

  # If evertyhing is already a single gene, return the original eset
  if (!any(grepl(separator, genes_to_map_raw))) {
    return(base_eset)
  }

  single_genes <- sapply(
    strsplit(genes_to_map_raw, separator), function(x) trimws(x[1])
  )
  fData(base_eset)[[gene_col]] <- single_genes
  return(base_eset)
}

#' Removes genes with empty gene IDs from an ExpressionSet object in a given
#' column.
#'
#' The function takes an ExpressionSet object (`eset`) and the name
#' of the column containing the gene IDs (`gene_col`) and removes any genes
#' with empty gene IDs.
#'
#' @param eset An ExpressionSet object.
#' @param gene_col The name of the column in the featureData of the
#' ExpressionSet object that contains the gene IDs.
#'
#' @return An ExpressionSet object with the empty genes removed.
#'
#' @examples
#' eset <- remove_empty_genes(eset, "ENTREZ_GENE_ID")
#'
#' @export
remove_empty_genes <- function(eset, gene_col) {
  base_eset <- eset
  gene_ids <- fData(base_eset)[[gene_col]]

  # Replace the missing value with NA and remove them
  gene_ids[gene_ids == ""] <- NA
  no_gene_id <- which(is.na(gene_ids))

  return(base_eset[-no_gene_id, ])
}

#' Aggregate expression data for genes with same gene IDs
#'
#' Aggregates the expression data for genes with the same gene IDs
#' and returns a new ExpressionSet object with the aggregated data.
#'
#' @param eset An ExpressionSet object containing gene expression data.
#' @param gene_col A character string specifying the column name in the fData
#' of the ExpressionSet object that contains the Entrez IDs of the genes.
#' @param aggregate_fun A function to aggregate the expression data for genes.
#' Default is 'max'.
#'
#' @return A new ExpressionSet object with the aggregated expression data.
#'
#' @examples
#' # Aggregate expression data by Entrez IDs
#' new_eset <- aggregate_expression(eset, "EntrezID")
#'
#' @importFrom Biobase ExpressionSet
#' @importFrom dplyr bind_cols group_by summarise
#' @importFrom dplyr ungroup across where
#' @importFrom rlang syms
aggregate_expression <- function(eset, gene_col, aggregate_fun = max) {
  base_eset <- eset

  # Create a temporary data frame combining fData and exprs
  tempDF <- bind_cols(fData(base_eset), as.data.frame(exprs(base_eset)))

  # Aggregate expression data for genes with same Entrez IDs
  aggregated_data_all <- tempDF %>%
    group_by(!!!syms(gene_col)) %>%
    summarise(across(where(is.numeric), aggregate_fun))

  aggregated_data <- aggregated_data_all[, c(gene_col, sampleNames(base_eset))]

  # Get the feature data for the remaining genes
  gene_idxs <- match(aggregated_data[[gene_col]], fData(base_eset)[[gene_col]])
  filtered_feature_data <- featureData(base_eset)[gene_idxs, ]

  new_eset <- ExpressionSet(
    assayData = as.matrix(aggregated_data[, -1]),
    phenoData = AnnotatedDataFrame(pData(base_eset)),
    featureData = filtered_feature_data
  )

  return(new_eset)
}


map_to_ensembl <- function(eset, gene_col, attribute, mart, handle = NULL) {
  base_eset <- eset
  gene_ids <- fData(base_eset)[[gene_col]]


  mapping <- biomaRt::getBM(
    attributes = c(attribute, "ensembl_gene_id"),
    values = gene_ids,
    filters = attribute,
    mart = mart
  )

  # Get index of genes with no mapping
  no_mapping <- which(!gene_ids %in% mapping[[attribute]])

  # Remove genes with no mapping
  mapped_genes_eset <- base_eset[-no_mapping, ]

  # Check if there are any genes left
  if (nrow(mapped_genes_eset) == 0) {
    print("No genes were mapped to Ensembl IDs")
  }

  unique_mapping <- mapping[!duplicated(mapping[[attribute]]), ]

  fData(mapped_genes_eset)[["ensembl_id"]] <- sapply(
    fData(mapped_genes_eset)[[gene_col]], function(x) {
      mapping$ensembl_gene_id[which(unique_mapping[[attribute]] == x)]
    }
  )


  # Duplicate expression if same gene map to different Ensembl ID
  duplicated_genes <- mapping[duplicated(mapping[[attribute]]), attribute]
  duplicated_genes_index <-
    which(featureNames(mapped_genes_eset) %in% duplicated_genes)

  if (length(duplicated_genes_index) == 0) {
    return(mapped_genes_eset)
  }

  duplicated_genes_exprs <- exprs(mapped_genes_eset)[duplicated_genes_index, ]

  # rename the rows with Ensembl ID
  rownames(duplicated_genes_exprs) <-
    mapping$ensembl_gene_id[duplicated_genes_index]

  new_eset <- ExpressionSet(
    assayData = rbind(exprs(mapped_genes_eset), duplicated_genes_exprs),
    phenoData = pData(mapped_genes_eset),
    featureData = fData(mapped_genes_eset)
  )

  return(new_eset)
}
