@testsnippet TestingEset begin
    using ExpressionData

    eset = load_eset(joinpath(@__DIR__, "data/test_eset.rds"))
end

@testitem "valid eset" setup = [TestingEset] begin
    annotated_eset = map_to_ensembl(eset, "entrezgene_id"; gene_col="ENTREZ_GENE_ID")

    # Test if the mapping was successful
    feature_names(annotated_eset) == feature_names(eset)

    @test "ensembl_id" in keys(feature_names(annotated_eset))
end
