push!(LOAD_PATH,"../src/")
push!(LOAD_PATH,"../../")

using Documenter, AComVarDesign

makedocs(sitename="My Documentation")

deploydocs(
    repo = "github.com/JoshuaLukemire/AComVarDesign.git",
)
