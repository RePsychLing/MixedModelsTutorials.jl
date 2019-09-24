using Documenter, Tutorials

makedocs(;
    modules=[Tutorials],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/bates/Tutorials.jl/blob/{commit}{path}#L{line}",
    sitename="Tutorials.jl",
    authors="Douglas Bates",
    assets=String[],
)

deploydocs(;
    repo="github.com/bates/Tutorials.jl",
)
