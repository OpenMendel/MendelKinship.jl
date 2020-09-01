using Documenter
using MendelKinship

makedocs(
    sitename = "MendelKinship",
    format = Documenter.HTML(),
    modules = [MendelKinship],
    pages = [
        "Home" => "index.md",
        "Introduction" => "man/introduction.md",
        "Tutorial" => "man/KinshipTutorial.md",
        "Keyword Options" => "man/keywords.md",
        "API" => "man/api.md"
    ]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo   = "github.com/OpenMendel/MendelKinship.jl.git",
    target = "build"
)
