
using Tensors
import Gmsh: gmsh

const lobatto3 = ([-1.0,0.0,0.0,
                    0.0,0.0,0.0,
                    1.0,0.0,0.0],[1/3,4/3,1/3])

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
    ∂²𝗠∂x² = (0,zeros(21))
    ∂²𝗠∂y² = (0,zeros(21))
    ∂²𝗠∂x∂y = (0,zeros(21))
    ∂³𝗠∂x³ = (0,zeros(21))
    ∂³𝗠∂x²∂y = (0,zeros(21))
    ∂³𝗠∂x∂y² = (0,zeros(21))
    ∂³𝗠∂y³ = (0,zeros(21))
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠, :∂²𝝭∂x²=>:𝑠, :∂²𝝭∂x∂y=>:𝑠, :∂²𝝭∂y²=>:𝑠)
    push!(elements["Ω"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂²𝗠∂x²=>∂²𝗠∂x², :∂²𝗠∂y²=>∂²𝗠∂y², :∂²𝗠∂x∂y=>∂²𝗠∂x∂y)
    push!(elements["Γᵇ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠, :∂²𝝭∂x²=>:𝑠, :∂²𝝭∂x∂y=>:𝑠, :∂²𝝭∂y²=>:𝑠, :∂³𝝭∂x³=>:𝑠, :∂³𝝭∂x²∂y=>:𝑠, :∂³𝝭∂x∂y²=>:𝑠, :∂³𝝭∂y³=>:𝑠)
    push!(elements["Γᵇ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂²𝗠∂x²=>∂²𝗠∂x², :∂²𝗠∂y²=>∂²𝗠∂y², :∂²𝗠∂x∂y=>∂²𝗠∂x∂y, :∂³𝗠∂x³=>∂³𝗠∂x³, :∂³𝗠∂x²∂y=>∂³𝗠∂x²∂y, :∂³𝗠∂x∂y²=>∂³𝗠∂x∂y², :∂³𝗠∂y³=>∂³𝗠∂y³)
    push!(elements["Γʳ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γʳ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γᵗ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠, :∂²𝝭∂x²=>:𝑠, :∂²𝝭∂x∂y=>:𝑠, :∂²𝝭∂y²=>:𝑠, :∂³𝝭∂x³=>:𝑠, :∂³𝝭∂x²∂y=>:𝑠, :∂³𝝭∂x∂y²=>:𝑠, :∂³𝝭∂y³=>:𝑠)
    push!(elements["Γᵗ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂²𝗠∂x²=>∂²𝗠∂x², :∂²𝗠∂y²=>∂²𝗠∂y², :∂²𝗠∂x∂y=>∂²𝗠∂x∂y, :∂³𝗠∂x³=>∂³𝗠∂x³, :∂³𝗠∂x²∂y=>∂³𝗠∂x²∂y, :∂³𝗠∂x∂y²=>∂³𝗠∂x∂y², :∂³𝗠∂y³=>∂³𝗠∂y³)
    push!(elements["Γˡ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠, :∂²𝝭∂x²=>:𝑠, :∂²𝝭∂x∂y=>:𝑠, :∂²𝝭∂y²=>:𝑠, :∂³𝝭∂x³=>:𝑠, :∂³𝝭∂x²∂y=>:𝑠, :∂³𝝭∂x∂y²=>:𝑠, :∂³𝝭∂y³=>:𝑠)
    push!(elements["Γˡ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂²𝗠∂x²=>∂²𝗠∂x², :∂²𝗠∂y²=>∂²𝗠∂y², :∂²𝗠∂x∂y=>∂²𝗠∂x∂y, :∂³𝗠∂x³=>∂³𝗠∂x³, :∂³𝗠∂x²∂y=>∂³𝗠∂x²∂y, :∂³𝗠∂x∂y²=>∂³𝗠∂x∂y², :∂³𝗠∂y³=>∂³𝗠∂y³)
    push!(elements["𝐴"], :𝝭=>:𝑠)
    push!(elements["𝐴"], :𝗠=>𝗠)

    # gmsh.finalize()
    return elements, nodes
end

function import_roof_mix(filename::String,n)
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

    integrationScheme = ([0.5,0.5,0.0,
                          0.5,0.0,0.5,
                          0.0,0.5,0.5],[1/3,1/3,1/3])
    # elements["Ωₚ"] = getMacroElementsForTriangles(entities["Ω"],PiecewisePolynomial{:Linear2D},integrationScheme,n)
    # elements["Ω"] = getElements(nodes, entities["Ω"], type, integrationScheme, sp)

    integrationScheme = lobatto7
    elements["Γₚ"] = getMacroBoundaryElementsForTriangles(entities["Γ"],entities["Ω"],PiecewisePolynomial{:Linear2D},integrationScheme,n)
    elements["Γ"] = getElements(nodes, entities["Γ"], type, integrationScheme, sp, normal = true)
    elements["Γᵇ"] = getElements(nodes, entities["Γᵇ"], type, integrationScheme, sp, normal = true)
    elements["Γʳ"] = getElements(nodes, entities["Γʳ"], type, integrationScheme, sp, normal = true)
    elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"], type, integrationScheme, sp, normal = true)
    elements["Γˡ"] = getElements(nodes, entities["Γˡ"], type, integrationScheme, sp, normal = true)

    elements["Ωₚ"] = getMacroElementsForTriangles(entities["Ω"],PiecewisePolynomial{:Linear2D},integrationOrder,n)
    elements["Ω"] = getElements(nodes, entities["Ω"], type, integrationOrder, sp)
    # elements["Γₚ"] = getMacroBoundaryElementsForTriangles(entities["Γ"],entities["Ω"],PiecewisePolynomial{:Linear2D},integrationOrder,n)
    # elements["Γ"] = getElements(nodes, entities["Γ"], type, integrationOrder, sp, normal = true)
    # elements["Γᵇ"] = getElements(nodes, entities["Γᵇ"], type, integrationOrder, sp, normal = true)
    # elements["Γʳ"] = getElements(nodes, entities["Γʳ"], type, integrationOrder, sp, normal = true)
    # elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"], type, integrationOrder, sp, normal = true)
    # elements["Γˡ"] = getElements(nodes, entities["Γˡ"], type, integrationOrder, sp, normal = true)

    elements["𝐴"] = getElements(nodes, entities["𝐴"], type, integrationOrder, sp)
    𝗠 = (0,zeros(21))
    ∂𝗠∂x = (0,zeros(21))
    ∂𝗠∂y = (0,zeros(21))
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ω"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Ωₚ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠, :∂²𝝭∂x²=>:𝑠, :∂²𝝭∂x∂y=>:𝑠, :∂²𝝭∂y²=>:𝑠)
    push!(elements["Γ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γₚ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
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

    # gmsh.finalize()
    return elements, nodes
end

prescribleBoundary = quote
    prescribe!(elements["Ω"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ω"],:a²²=>(ξ¹,ξ²,ξ³)->cs.a²²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ω"],:a¹²=>(ξ¹,ξ²,ξ³)->cs.a¹²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ω"],:𝒂₁₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Ω"],:𝒂₁₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Ω"],:𝒂₁₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Ω"],:𝒂₂₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Ω"],:𝒂₂₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Ω"],:𝒂₂₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Ω"],:𝒂₃₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Ω"],:𝒂₃₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Ω"],:𝒂₃₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Ω"],:Γ¹₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ω"],:Γ¹₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ω"],:Γ¹₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ω"],:Γ²₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ω"],:Γ²₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ω"],:Γ²₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ω"],:b₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Ω"],:b₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Ω"],:b₃=>(ξ¹,ξ²,ξ³)->b₃)
    prescribe!(elements["Γᵇ"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:a²²=>(ξ¹,ξ²,ξ³)->cs.a²²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:a¹²=>(ξ¹,ξ²,ξ³)->cs.a¹²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵇ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵇ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵇ"],:n₁₁=>(ξ¹,ξ²,ξ³)->1.0)
    prescribe!(elements["Γᵇ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵇ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵇ"],:n₂₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵇ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵇ"],:n₃₃=>(ξ¹,ξ²,ξ³)->1.0)
    prescribe!(elements["Γᵇ"],:𝒂₁₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γᵇ"],:𝒂₁₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γᵇ"],:𝒂₁₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γᵇ"],:𝒂₂₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γᵇ"],:𝒂₂₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γᵇ"],:𝒂₂₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γᵇ"],:𝒂₃₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γᵇ"],:𝒂₃₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γᵇ"],:𝒂₃₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γᵗ"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:a²²=>(ξ¹,ξ²,ξ³)->cs.a²²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:a¹²=>(ξ¹,ξ²,ξ³)->cs.a¹²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵗ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵗ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵗ"],:n₁₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵗ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵗ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵗ"],:n₂₂=>(ξ¹,ξ²,ξ³)->1.0)
    prescribe!(elements["Γᵗ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵗ"],:n₃₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵗ"],:𝒂₃₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γᵗ"],:𝒂₃₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γᵗ"],:𝒂₃₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γᵗ"],:θ=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γᵗ"],:𝒂₁₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γᵗ"],:𝒂₁₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γᵗ"],:𝒂₁₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γᵗ"],:𝒂₂₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γᵗ"],:𝒂₂₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γᵗ"],:𝒂₂₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γᵗ"],:𝒂₃₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γᵗ"],:𝒂₃₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γᵗ"],:𝒂₃₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γˡ"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:a²²=>(ξ¹,ξ²,ξ³)->cs.a²²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:a¹²=>(ξ¹,ξ²,ξ³)->cs.a¹²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₁₁=>(ξ¹,ξ²,ξ³)->1.0)
    prescribe!(elements["Γˡ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₂₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₃₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:θ=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:𝒂₃₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γˡ"],:𝒂₃₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γˡ"],:𝒂₃₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γˡ"],:𝒂₁₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γˡ"],:𝒂₁₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γˡ"],:𝒂₁₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γˡ"],:𝒂₂₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γˡ"],:𝒂₂₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γˡ"],:𝒂₂₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γˡ"],:𝒂₃₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γˡ"],:𝒂₃₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γˡ"],:𝒂₃₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
end

prescribleForMix = quote
    prescribe!(elements["Ωₚ"],:a₁₁=>(ξ¹,ξ²,ξ³)->cs.a₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:a₂₂=>(ξ¹,ξ²,ξ³)->cs.a₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:a₁₂=>(ξ¹,ξ²,ξ³)->cs.a₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:𝒂¹₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Ωₚ"],:𝒂¹₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Ωₚ"],:𝒂¹₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Ωₚ"],:𝒂²₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Ωₚ"],:𝒂²₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Ωₚ"],:𝒂²₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Ωₚ"],:Γ¹₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ¹₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ¹₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ²₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ²₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ²₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ¹₁₁₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ¹₁₁₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ¹₁₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ¹₁₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ¹₂₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ¹₂₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ²₁₁₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ²₁₁₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ²₁₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ²₁₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ²₂₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:Γ²₂₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    
    prescribe!(elements["Γₚ"],:Γ¹₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ¹₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ¹₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ²₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ²₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ²₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:a²²=>(ξ¹,ξ²,ξ³)->cs.a²²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:a¹²=>(ξ¹,ξ²,ξ³)->cs.a¹²(Vec{3}((ξ¹,ξ²,ξ³))))
end