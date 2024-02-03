
using Tensors, BenchmarkExample
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

const trilobatto13 = ([1.0000000000000000,0.0000000000000000,0.0,
                       0.0000000000000000,1.0000000000000000,0.0,
                       0.0000000000000000,0.0000000000000000,0.0,
                       0.0000000000000000,0.5000000000000000,0.0,
                       0.5000000000000000,0.0000000000000000,0.0,
                       0.5000000000000000,0.5000000000000000,0.0,
                       0.0000000000000000,0.8273268353539885,0.0,
                       0.0000000000000000,0.1726731646460114,0.0,
                       0.1726731646460114,0.0000000000000000,0.0,
                       0.8273268353539885,0.0000000000000000,0.0,
                       0.8273268353539885,0.1726731646460114,0.0,
                       0.1726731646460114,0.8273268353539885,0.0,
                       0.3333333333333333,0.3333333333333333,0.0],
                       [-0.0277777777777778,
                        -0.0277777777777778,
                        -0.0277777777777778,
                         0.0296296296296297,
                         0.0296296296296297,
                         0.0296296296296297,
                         0.0907407407407407,
                         0.0907407407407407,
                         0.0907407407407407,
                         0.0907407407407407,
                         0.0907407407407407,
                         0.0907407407407407,
                         0.4500000000000000])

function import_roof_gauss(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
    integrationOrder = 3
    nₘ = 55
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    cs = BenchmarkExample.cylindricalCoordinate(BenchmarkExample.ScordelisLoRoof.𝑅)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getCurvedElements(nodes, entities["Ω"], type, cs, integrationOrder, sp)
    elements["Γᵇ"] = getCurvedElements(nodes, entities["Γᵇ"], type, cs, integrationOrder, sp)
    elements["Γʳ"] = getCurvedElements(nodes, entities["Γʳ"], type, cs, integrationOrder, sp)
    elements["Γᵗ"] = getCurvedElements(nodes, entities["Γᵗ"], type, cs, integrationOrder, sp)
    elements["Γˡ"] = getCurvedElements(nodes, entities["Γˡ"], type, cs, integrationOrder, sp)
    elements["𝐴"] = getElements(nodes, entities["𝐴"], type, integrationOrder, sp)
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    ∂²𝗠∂x² = (0,zeros(nₘ))
    ∂²𝗠∂y² = (0,zeros(nₘ))
    ∂²𝗠∂x∂y = (0,zeros(nₘ))
    ∂³𝗠∂x³ = (0,zeros(nₘ))
    ∂³𝗠∂x²∂y = (0,zeros(nₘ))
    ∂³𝗠∂x∂y² = (0,zeros(nₘ))
    ∂³𝗠∂y³ = (0,zeros(nₘ))
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

    gmsh.finalize()
    return elements, nodes
end

function import_roof_mix(filename::String,n)
    gmsh.initialize()
    gmsh.open(filename)
    integrationOrder = 3    
    # type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
    # nₘ = 21
    type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
    nₘ = 55
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    cs = BenchmarkExample.cylindricalCoordinate(BenchmarkExample.ScordelisLoRoof.𝑅)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    integrationScheme = trilobatto13
    elements["Ωₚ"] = getCurvedPiecewiseElements(entities["Ω"], PiecewisePolynomial{:Linear2D}, cs, integrationScheme)
    elements["Ω"] = getCurvedElements(nodes, entities["Ω"], type, cs, integrationScheme, sp)

    integrationScheme = lobatto5
    elements["Γₚ"] = getCurvedPiecewiseElements(entities["Γ"],PiecewisePolynomial{:Linear2D}, cs, integrationScheme,3)
    elements["Γ"] = getCurvedElements(nodes, entities["Γ"], type, cs, integrationScheme, sp)
    elements["Γᵇ"] = getCurvedElements(nodes, entities["Γᵇ"], type, cs, integrationScheme, sp)
    elements["Γʳ"] = getCurvedElements(nodes, entities["Γʳ"], type, cs, integrationScheme, sp)
    elements["Γᵗ"] = getCurvedElements(nodes, entities["Γᵗ"], type, cs, integrationScheme, sp)
    elements["Γˡ"] = getCurvedElements(nodes, entities["Γˡ"], type, cs, integrationScheme, sp)

    # elements["Ωₚ"] = getMacroElements(entities["Ω"],PiecewisePolynomial{:Linear2D},integrationOrder,n,nₕ=2,nₐ=2)
    # elements["Ω"] = getElements(nodes, entities["Ω"], type, integrationOrder, sp)
    # elements["Γₚ"] = getMacroBoundaryElementsForTriangles(entities["Γ"],entities["Ω"],PiecewisePolynomial{:Linear2D},integrationOrder,n)
    # elements["Γ"] = getElements(nodes, entities["Γ"], type, integrationOrder, sp, normal = true)
    # elements["Γᵇ"] = getElements(nodes, entities["Γᵇ"], type, integrationOrder, sp, normal = true)
    # elements["Γʳ"] = getElements(nodes, entities["Γʳ"], type, integrationOrder, sp, normal = true)
    # elements["Γᵗ"] = getElements(nodes, entities["Γᵗ"], type, integrationOrder, sp, normal = true)
    # elements["Γˡ"] = getElements(nodes, entities["Γˡ"], type, integrationOrder, sp, normal = true)

    elements["𝐴"] = getElements(nodes, entities["𝐴"], type, integrationOrder, sp)
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    ∂²𝗠∂x² = (0,zeros(nₘ))
    ∂²𝗠∂x∂y = (0,zeros(nₘ))
    ∂²𝗠∂y² = (0,zeros(nₘ))
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠, :∂²𝝭∂x²=>:𝑠, :∂²𝝭∂x∂y=>:𝑠, :∂²𝝭∂y²=>:𝑠)
    push!(elements["Ω"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂²𝗠∂x²=>∂²𝗠∂x², :∂²𝗠∂x∂y=>∂²𝗠∂x∂y, :∂²𝗠∂y²=>∂²𝗠∂y²)
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
    prescribe!(elements["Γᵇ"],:b₁₁=>(ξ¹,ξ²,ξ³)->cs.b₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:b₂₂=>(ξ¹,ξ²,ξ³)->cs.b₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:b₁₂=>(ξ¹,ξ²,ξ³)->cs.b₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ¹₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ¹₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ¹₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ²₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ²₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ²₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ¹₁₁₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ¹₁₁₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ¹₁₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ¹₁₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ¹₂₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ¹₂₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ²₁₁₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ²₁₁₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ²₁₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ²₁₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ²₂₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵇ"],:Γ²₂₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
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
    prescribe!(elements["Γᵇ"],:𝒂¹₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γᵇ"],:𝒂¹₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γᵇ"],:𝒂¹₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γᵇ"],:𝒂²₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γᵇ"],:𝒂²₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γᵇ"],:𝒂²₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γᵇ"],:𝒂³₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂³(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γᵇ"],:𝒂³₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂³(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γᵇ"],:𝒂³₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂³(Vec{3}((ξ¹,ξ²,ξ³)))[3])

    prescribe!(elements["Γᵗ"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:a²²=>(ξ¹,ξ²,ξ³)->cs.a²²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:a¹²=>(ξ¹,ξ²,ξ³)->cs.a¹²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:b₁₁=>(ξ¹,ξ²,ξ³)->cs.b₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:b₂₂=>(ξ¹,ξ²,ξ³)->cs.b₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:b₁₂=>(ξ¹,ξ²,ξ³)->cs.b₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ¹₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ¹₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ¹₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ²₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ²₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ²₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ¹₁₁₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ¹₁₁₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ¹₁₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ¹₁₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ¹₂₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ¹₂₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ²₁₁₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ²₁₁₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ²₁₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ²₁₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ²₂₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γᵗ"],:Γ²₂₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
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
    prescribe!(elements["Γᵗ"],:𝒂¹₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γᵗ"],:𝒂¹₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γᵗ"],:𝒂¹₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γᵗ"],:𝒂²₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γᵗ"],:𝒂²₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γᵗ"],:𝒂²₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γᵗ"],:𝒂³₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂³(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γᵗ"],:𝒂³₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂³(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γᵗ"],:𝒂³₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂³(Vec{3}((ξ¹,ξ²,ξ³)))[3])

    prescribe!(elements["Γˡ"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:a²²=>(ξ¹,ξ²,ξ³)->cs.a²²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:a¹²=>(ξ¹,ξ²,ξ³)->cs.a¹²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:b₁₁=>(ξ¹,ξ²,ξ³)->cs.b₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:b₂₂=>(ξ¹,ξ²,ξ³)->cs.b₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:b₁₂=>(ξ¹,ξ²,ξ³)->cs.b₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ¹₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ¹₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ¹₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ²₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ²₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ²₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ¹₁₁₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ¹₁₁₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ¹₁₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ¹₁₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ¹₂₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ¹₂₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ²₁₁₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ²₁₁₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ²₁₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ²₁₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ²₂₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:Γ²₂₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
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
    prescribe!(elements["Γˡ"],:𝒂¹₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γˡ"],:𝒂¹₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γˡ"],:𝒂¹₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γˡ"],:𝒂²₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γˡ"],:𝒂²₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γˡ"],:𝒂²₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γˡ"],:𝒂³₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂³(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γˡ"],:𝒂³₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂³(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γˡ"],:𝒂³₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂³(Vec{3}((ξ¹,ξ²,ξ³)))[3])
end

prescribleForMix = quote
    prescribe!(elements["Ωₚ"],:a₁₁=>(ξ¹,ξ²,ξ³)->cs.a₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:a₂₂=>(ξ¹,ξ²,ξ³)->cs.a₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:a₁₂=>(ξ¹,ξ²,ξ³)->cs.a₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:b₁₁=>(ξ¹,ξ²,ξ³)->cs.b₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:b₂₂=>(ξ¹,ξ²,ξ³)->cs.b₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:b₁₂=>(ξ¹,ξ²,ξ³)->cs.b₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Ωₚ"],:𝒂₁₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Ωₚ"],:𝒂₁₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Ωₚ"],:𝒂₁₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Ωₚ"],:𝒂₂₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Ωₚ"],:𝒂₂₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Ωₚ"],:𝒂₂₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Ωₚ"],:𝒂₃₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Ωₚ"],:𝒂₃₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Ωₚ"],:𝒂₃₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Ωₚ"],:𝒂¹₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Ωₚ"],:𝒂¹₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Ωₚ"],:𝒂¹₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Ωₚ"],:𝒂²₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Ωₚ"],:𝒂²₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Ωₚ"],:𝒂²₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[3])
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
    prescribe!(elements["Ωₚ"],:𝒂³₁₁₍₁₎=>(ξ¹,ξ²,ξ³)->cs.∂₁₁𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Ωₚ"],:𝒂³₁₁₍₂₎=>(ξ¹,ξ²,ξ³)->cs.∂₁₁𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Ωₚ"],:𝒂³₁₁₍₃₎=>(ξ¹,ξ²,ξ³)->cs.∂₁₁𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Ωₚ"],:𝒂³₂₂₍₁₎=>(ξ¹,ξ²,ξ³)->cs.∂₂₂𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Ωₚ"],:𝒂³₂₂₍₂₎=>(ξ¹,ξ²,ξ³)->cs.∂₂₂𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Ωₚ"],:𝒂³₂₂₍₃₎=>(ξ¹,ξ²,ξ³)->cs.∂₂₂𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Ωₚ"],:𝒂³₁₂₍₁₎=>(ξ¹,ξ²,ξ³)->cs.∂₁₂𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Ωₚ"],:𝒂³₁₂₍₂₎=>(ξ¹,ξ²,ξ³)->cs.∂₁₂𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Ωₚ"],:𝒂³₁₂₍₃₎=>(ξ¹,ξ²,ξ³)->cs.∂₁₂𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    
    prescribe!(elements["Γₚ"],:Γ¹₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ¹₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ¹₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ²₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ²₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ²₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:a²²=>(ξ¹,ξ²,ξ³)->cs.a²²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:a¹²=>(ξ¹,ξ²,ξ³)->cs.a¹²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:b₁₁=>(ξ¹,ξ²,ξ³)->cs.b₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:b₂₂=>(ξ¹,ξ²,ξ³)->cs.b₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:b₁₂=>(ξ¹,ξ²,ξ³)->cs.b₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:𝒂₁₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γₚ"],:𝒂₁₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γₚ"],:𝒂₁₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γₚ"],:𝒂₂₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γₚ"],:𝒂₂₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γₚ"],:𝒂₂₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γₚ"],:𝒂₃₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γₚ"],:𝒂₃₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γₚ"],:𝒂₃₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γₚ"],:𝒂¹₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γₚ"],:𝒂¹₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γₚ"],:𝒂¹₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γₚ"],:𝒂²₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γₚ"],:𝒂²₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γₚ"],:𝒂²₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[3])
end
