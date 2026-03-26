using Documenter
using SignatureTensors 

makedocs(
    sitename = "SignatureTensors.jl",
    modules = [SignatureTensors],
    pages = [
        "Home" => "index.md",
        "Documentation" => "api.md",
    ],
    checkdocs = :warn,  
)



deploydocs(
    repo = "github.com/leonardSchmitz/signature-tensors-in-OSCAR.git",
)

