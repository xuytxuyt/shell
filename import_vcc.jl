
import Gmsh: gmsh

const lobatto3 = ([-1.0,0.0,0.0,
                    0.0,0.0,0.0,
                    1.0,0.0,0.0],[1/3,4/3,1/3])

const lobatto5 = ([-1.0,0.0,0.0,
                   -(3/7)^0.5,0.0,0.0,
                    0.0,0.0,0.0,
                    (3/7)^0.5,0.0,0.0,
                    1.0,0.0,0.0],[1/10,49/90,32/45,49/90,1/10])

const lobatto7 = ([-1.0,0.0,0.0,
                   -(5/11+2/11*(5/3)^0.5)^0.5,0.0,0.0,
                   -(5/11-2/11*(5/3)^0.5)^0.5,0.0,0.0,
                    0.0,0.0,0.0,
                    (5/11-2/11*(5/3)^0.5)^0.5,0.0,0.0,
                    (5/11+2/11*(5/3)^0.5)^0.5,0.0,0.0,
                    1.0,0.0,0.0],
                  [1/21,
                   (124-7*15^0.5)/350,
                   (124+7*15^0.5)/350,
                   256/525,
                   (124+7*15^0.5)/350,
                   (124-7*15^0.5)/350,
                   1/21])

function import_vcc(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    cs = BenchmarkExample.cylindricalCoordinate(25.0)
    # cs = BenchmarkExample.sphericalCoordinate(10.0)
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
    # elements["Γ"] = getCurvedElements(nodes, entities["Γ"], cs, integrationOrder)
    elements["Γ"] = getCurvedElements(nodes, entities["Γ"], cs, lobatto5)

    return elements, nodes
end