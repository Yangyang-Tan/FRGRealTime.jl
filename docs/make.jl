using FRGRealTime
using Documenter

DocMeta.setdocmeta!(FRGRealTime, :DocTestSetup, :(using FRGRealTime); recursive=true)

makedocs(;
    modules=[FRGRealTime],
    authors="Yang-yang Tan, Yong-rui Chen, Wei-jie Fu",
    repo="https://github.com/Yangyang-Tan/FRGRealTime.jl/blob/{commit}{path}#{line}",
    sitename="FRGRealTime.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Yangyang-Tan.github.io/FRGRealTime.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Example" => "markdown/Example.md"
    ],
)

deploydocs(;
    repo="github.com/Yangyang-Tan/FRGRealTime.jl",
)
