using Documenter
using SignatureTensors
using DocumenterCitations

# Bibliografía

bib = CitationBibliography(joinpath(@__DIR__, "refs.bib"), style = :numeric)
# Generar la documentación
makedocs(
    sitename = "SignatureTensors.jl",
    modules = [SignatureTensors],
    pages = [
        "Home" => "index.md",
        "Documentation" => "api.md",
        "References" => "references.md",
    ],

    checkdocs = :warn,
    plugins = [bib],
    highlightsig = true,

    meta = Dict(
        :CollapsedDocStrings => true
    ),
    
    format = Documenter.HTML(
        prettyurls = false
    )
)

build_dir = joinpath(@__DIR__, "build")

for (root, dirs, files) in walkdir(build_dir)
    for file in files
        if endswith(file, ".html")
            path = joinpath(root, file)
            html = read(path, String)

            html = replace( html, r"(<summary[^>]*>.*?<code>)SignatureTensors\.(.*?</code>.*?</summary>)"s  => s"\1\2" )

            write(path, html)
        end
    end
end


deploydocs(
    repo = "github.com/leonardSchmitz/SignatureTensors.jl",
    branch = "gh-pages",
    devurl = "docs",
    forcepush = true
)

