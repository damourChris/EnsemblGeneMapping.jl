module EnsemblGeneMapping

using Comonicon
using RCall
using ExpressionData

@main function map_to_ensembl(eset::String, attribute::String;
                              output_file::Union{Nothing,String}=nothing,
                              config_file::Union{Nothing,String}=nothing)
    # Parse the arguments into a config struct
    config::Config = load_config(config_file)

    output_file = isnothing(output_file) ? eset : output_file

    eset = load_eset(eset)

    (; gene_col, mart_id, mart_dataset) = config

    # Get the R utils files path
    r_utils_path = joinpath(@__DIR__, "utils.r")

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

    return save_eset(eset, output_file)
end

end
