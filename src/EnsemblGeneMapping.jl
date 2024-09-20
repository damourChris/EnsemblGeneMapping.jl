module EnsemblGeneMapping

using ExpressionData
using RCall

export map_to_ensembl

# Get the R utils files path
const r_utils_path = joinpath(@__DIR__, "utils.r")

function map_to_ensembl(eset::ExpressionSet, attribute::String; 
        gene_col::String = "gene_symbol",  
        mart_id::String = "ensembl", 
        mart_dataset::String = "hsapiens_gene_ensembl",
    )

    @rput eset gene_col mart_id mart_dataset attribute r_utils_path

    R"""
    source(r_utils_path)
    mart <- biomaRt::useMart(mart_id, dataset = mart_dataset)

    annotated_eset <-
        remove_empty_genes(eset, gene_col = gene_col) %>%
        extract_first_gene_symbol(gene_col = gene_col) %>%
        aggregate_expression(gene_col = gene_col) %>%
        map_to_ensembl(
            gene_col = gene_col,
            attribute = attribute,
            mart = mart
            ) %>%
        aggregate_expression(gene_col = "ensembl_id")
    """

    eset = @rget annotated_eset

    return eset
end

end
