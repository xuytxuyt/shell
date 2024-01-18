
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
    𝗠 = (0,zeros(21))
    ∂𝗠∂x = (0,zeros(21))
    ∂𝗠∂y = (0,zeros(21))
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠, :∂²𝝭∂x²=>:𝑠, :∂²𝝭∂x∂y=>:𝑠, :∂²𝝭∂y²=>:𝑠,)
    push!(elements["Ω"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂²𝗠∂x²=>(0,zeros(21)), :∂²𝗠∂y²=>(0,zeros(21)), :∂²𝗠∂x∂y=>(0,zeros(21)),)
    push!(elements["Γᵇ"], :𝝭=>:𝑠)
    push!(elements["Γᵇ"], :𝗠=>𝗠)
    push!(elements["Γʳ"], :𝝭=>:𝑠)
    push!(elements["Γʳ"], :𝗠=>𝗠)
    push!(elements["Γᵗ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γᵗ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γˡ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γˡ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["𝐴"], :𝝭=>:𝑠)
    push!(elements["𝐴"], :𝗠=>𝗠)

    gmsh.finalize()
    return elements, nodes
end

function import_roof_mix(filename::String,n)
    gmsh.initialize()
    gmsh.open(filename)
    
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()

    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    integrationScheme = ([0.5,0.5,0.0,
                          0.5,0.0,0.5,
                          0.0,0.5,0.5],[1/3,1/3,1/3])
    elements["Ωₚ"] = getMacroElementsForTriangles(entities["Ω"],PiecewisePolynomial{:Linear2D},integrationScheme,n)

    integrationScheme = ([-1.0,0.0,0.0,
                           0.0,0.0,0.0,
                           1.0,0.0,0.0],[1/3,4/3,1/3])
    elements["Γₚ"] = getMacroBoundaryElementsForTriangles(entities["Γ"],entities["Ω"],PiecewisePolynomial{:Linear2D},integrationScheme,n)

    gmsh.finalize()
    return elements, nodes
end