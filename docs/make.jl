using Documenter, MendelKinship

ENV["DOCUMENTER_DEBUG"] = "true"

makedocs(
    format = :html,
    sitename = "MendelKinship",
    modules = [MendelKinship],
    authors = "Jeanette Papp",
    clean = true,
    debug = true,
    pages = [
        "Home" => "index.md",
        "API" => "man/api.md"
        "Tutorial" => "man/KinshipTutorial.md"
    ]
)

deploydocs(
    repo   = "github.com/OpenMendel/MendelKinship.jl.git",
    target = "build",
    julia  = "1.0",
    deps   = nothing,
    make   = nothing
)
