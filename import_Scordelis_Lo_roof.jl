
import Gmsh: gmsh

function import_roof_gauss(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"])
    elements["Γᵇ"] = getElements(nodes, entities["Γᵇ"])
    elements["Γʳ"] = getElements(nodes, entities["Γʳ"])
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"])
    elements["Γˡ"] = getElements(nodes, entities["Γˡ"])
    elements["𝐴"] = getElements(nodes, entities["𝐴"])

    gmsh.finalize()
    return elements, nodes
end