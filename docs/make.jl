using Documenter
using HydroForge

DocMeta.setdocmeta!(HydroForge, :DocTestSetup, :(using HydroForge); recursive=true)

makedocs(
    sitename = "HydroForge",
    modules = [HydroForge],
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        canonical = "https://hydroforge.github.io/HydroForge.jl",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
        "Getting Started" => [
            "Installation" => "installation.md",
            "Quick Start" => "quickstart.md",
        ],
        "User Guide" => [
            "Concepts" => "concepts.md",
            "Input Data" => "input_data.md",
            "Running Simulations" => "running.md",
            "Output Data" => "output_data.md",
        ],
        "API Reference" => [
            "Types" => "api/types.md",
            "Numerics" => "api/numerics.md",
            "Physics" => "api/physics.md",
            "IO" => "api/io.md",
        ],
        "Developer Guide" => [
            "Architecture" => "dev/architecture.md",
            "Contributing" => "dev/contributing.md",
        ],
    ],
    checkdocs = :exports,
    warnonly = [:missing_docs, :cross_references],
)

deploydocs(
    repo = "github.com/hydroforge/HydroForge.jl.git",
    devbranch = "main",
    push_preview = true,
)
