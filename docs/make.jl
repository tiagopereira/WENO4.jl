using WENO4
using Documenter

DocMeta.setdocmeta!(WENO4, :DocTestSetup, :(using WENO4); recursive=true)

makedocs(;
    modules=[WENO4],
    authors="Tiago M. D. Pereira",
    repo="https://github.com/tiagopereira/WENO4.jl/blob/{commit}{path}#{line}",
    sitename="WENO4.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://tiagopereira.github.io/WENO4.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/tiagopereira/WENO4.jl",
    devbranch="main",
)
