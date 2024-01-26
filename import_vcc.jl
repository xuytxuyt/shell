
import Gmsh: gmsh

function import_vcc(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    cs = BenchmarkExample.sphericalCoordinate(25.0)
    𝐽 = cs.𝐽
    𝒂₁ = cs.𝒂₁
    𝒂₂ = cs.𝒂₂
    𝒂₃ = cs.𝒂₃
    integrationOrder = 10
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getCurvedElements(nodes, entities["Ω"], integrationOrder,𝐽=𝐽)
    elements["Γ"] = getCurvedElements(nodes, entities["Γ"], integrationOrder, normal = true,𝒂₁=𝒂₁,𝒂₂=𝒂₂,𝒂₃=𝒂₃)

    return elements, nodes
end