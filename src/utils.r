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


  # Mapping assumptions
  # Diff OG id -> Same Ensembl ID -> Aggregate expression
  # Same OG id -> Same Ensembl ID -> Keep expression
  # Same OG id -> Diff Ensembl ID -> Duplicate expression

  #  Workflow: make a subset of the original Eset for each of the above 
  # cases and then combine them
  
  # Before anything, lets construct a df that indicates what to do with each gene
  # 1. If the gene is not mapped, we will remove it
  # 2. If the gene is mapped, we will keep it
  #  - a. If the gene is mapped to multiple ensembl ids, we will duplicate the expression 
  #  - b. If the gene is mapped to a single ensembl id, we will keep the expression as is
  #  - c. If multiple gene are mapped to the same ensembl id, we will aggregate the expression
  
  # Step 1: Identify unique mappings (one-to-one or many-to-one)
  unique_mappings <- mapping %>%
    group_by(!!sym(attribute)) %>%
    summarise(
      n = n(),
      ensembl_gene_id = list(unique(ensembl_gene_id))
    ) %>%
    filter(n == 1 | length(ensembl_gene_id) > 0)
  
  # Step 2: Process one-to-many mappings
  one_to_many <- unique_mappings %>%
    filter(length(ensembl_gene_id) == 1)
  
  if (nrow(one_to_many) > 0) {
    one_to_many_expanded <- one_to_many %>%
      unnest(ensembl_ids) %>%
      select(attribute, ensembl_id)
    
    df <- bind_rows(df, one_to_many_expanded)
  }
  
  # Step 3: Remove unmapped genes
  mapped_attributes <- unique(mapping[[attribute]])
    
    # Filter the original data frame to keep only mapped attributes
  new_fdata <- fData(base_eset) %>%
    filter(!!sym(gene_col) %in% mapped_attributes) 
  
  # Add the ensembl ids as a new column by matching the gene ids
  new_fdata$ensembl_id <- 
    mapping$ensembl_gene_id[match(new_fdata[[gene_col]], mapping[[attribute]])]
    
  # Step 4: Aggregate expressions for many-to-one mappings
    
  # First get the rows of the epxression matrix grouped by the ensembl id and remove the unmapped genes
  # which are stored in the new_fdata data frame
  gx_data <- exprs(base_eset)
  rownames(gx_data) <- fData(base_eset)[[gene_col]]
  gx_data <- gx_data[rownames(gx_data) %in% new_fdata[[gene_col]], ]

  # Aggregate the expression data for the genes with the same ensembl id
  aggregated_gx_data <- aggregate(
    gx_data, 
    by = list(ensembl_id = new_fdata$ensembl_id), 
    FUN = sum
  )
  
  # Remove the row names and the ensembl id column
  rownames(aggregated_gx_data) <- NULL
  aggregated_gx_data <- aggregated_gx_data[, -1]
  
  # Update the new_fdata data frame to remove the genes that were not mapped
  new_fdata <- new_fdata[match(rownames(aggregated_gx_data), new_fdata$ensembl_id), ]
  
  # Create a new ExpressionSet object with the aggregated expression data
  new_eset <- ExpressionSet(
    assayData = as.matrix(aggregated_gx_data),
    phenoData = AnnotatedDataFrame(pData(base_eset)),
    featureData = AnnotatedDataFrame(new_fdata)
  )
}