
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
    elements["Ω"] = getCurvedElements(nodes, entities["Ω"], integrationOrder,𝐽=𝐽)
    elements["Γ"] = getCurvedElements(nodes, entities["Γ"], integrationOrder,𝒂₁=𝒂₁,𝒂₂=𝒂₂,𝒂₃=𝒂₃,a₁₁=a₁₁,a₂₂=a₂₂,a₁₂=a₁₂,a¹¹=a¹¹,a²²=a²²,a¹²=a¹²)

    return elements, nodes
end