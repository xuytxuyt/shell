
import Gmsh: gmsh

function import_vcc(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    # cs = BenchmarkExample.cylindricalCoordinate(25.0)
    cs = BenchmarkExample.sphericalCoordinate(10.0)
    𝐽 = cs.𝐽
    𝒂₁ = cs.𝒂₁
    𝒂₂ = cs.𝒂₂
    𝒂₃ = cs.𝒂₃
    a₁₁ = cs.a₁₁
    a₂₂ = cs.a₂₂
    a₁₂ = cs.a₁₂
    a¹¹ = cs.a¹¹
    a²² = cs.a²²
    a¹² = cs.a¹²
    integrationOrder = 12
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getCurvedElements(nodes, entities["Ω"], cs, integrationOrder)
    elements["Γ"] = getCurvedElements(nodes, entities["Γ"], cs, integrationOrder)

    return elements, nodes
end