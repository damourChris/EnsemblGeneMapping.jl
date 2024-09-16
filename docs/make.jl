using EnsemblGeneMapping
using Documenter

DocMeta.setdocmeta!(EnsemblGeneMapping, :DocTestSetup, :(using EnsemblGeneMapping); recursive=true)

makedocs(;
    modules=[EnsemblGeneMapping],
    authors="Chris Damour",
    sitename="EnsemblGeneMapping.jl",
    format=Documenter.HTML(;
        canonical="https://damourChris.github.io/EnsemblGeneMapping.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/damourChris/EnsemblGeneMapping.jl",
    devbranch="main",
)
