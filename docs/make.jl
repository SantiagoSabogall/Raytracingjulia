using Documenter
using Raytracingjulia

makedocs(
    sitename = "Raytracingjulia",
    modules = [Raytracingjulia],
    pages = [
        "Inicio" => "index.md",
        "API" => "api.md",
    ],
    checkdocs = :none,
)
