@testsnippet TestingEset begin
    using ExpressionData

    eset = load_eset(joinpath(@__DIR__, "data/test_eset.rds"))
end

@testitem "valid eset" setup = [TestingEset] begin
    annotated_eset = map_to_ensembl(eset, "entrezgene_id"; gene_col="ENTREZ_GENE_ID")

    @test !ismissing(annotated_eset)
    @test !ismissing(feature_data(annotated_eset))
    @test !ismissing(expression_values(annotated_eset))
    @test !ismissing(phenotype_data(annotated_eset))
end
