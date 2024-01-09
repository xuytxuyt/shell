
import Gmsh: gmsh

function import_roof_gauss(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    integrationOrder = 8
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    sp = RegularGrid(x,y,z,n = 3,γ = 5)

    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getElements(nodes, entities["Ω"], type, integrationOrder, sp)
    elements["Γᵇ"] = getElements(nodes, entities["Γᵇ"], type, integrationOrder, sp, normal = true)
    elements["Γʳ"] = getElements(nodes, entities["Γʳ"], type, integrationOrder, sp, normal = true)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"], type, integrationOrder, sp, normal = true)
    elements["Γˡ"] = getElements(nodes, entities["Γˡ"], type, integrationOrder, sp, normal = true)
    elements["𝐴"] = getElements(nodes, entities["𝐴"], type, integrationOrder, sp)
    push!(elements["Ω"], :𝑠=>:∂²𝝭∂x², :𝑠=>:∂²𝝭∂x∂y, :𝑠=>:∂²𝝭∂y²)

    gmsh.finalize()
    return elements, nodes
end

function import_roof_mix(filename::String)
    gmsh.initialize()
    gmsh.open(filename)
    
    entities = getPhysicalGroups()
    addEdgeElements(entities["Ω"])

    nodes = get𝑿ᵢ()

    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    gmsh.finalize()
    return elements, nodes
end