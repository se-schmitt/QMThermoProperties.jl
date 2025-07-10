using QMIdealModels
using Documenter

DocMeta.setdocmeta!(QMIdealModels, :DocTestSetup, :(using QMIdealModels); recursive=true)

makedocs(;
    modules=[QMIdealModels],
    authors="Sebastian Schmitt <se.s.schmitt@gmail.com> and contributors",
    sitename="QMIdealModels.jl",
    format=Documenter.HTML(;
        canonical="https://se-schmitt.github.io/QMIdealModels.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/se-schmitt/QMIdealModels.jl",
    devbranch="main",
)
