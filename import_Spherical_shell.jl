
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
                   0.5*[-0.0277777777777778,
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

function import_spherical_gauss(filename::String)
    gmsh.initialize()
    gmsh.open(filename)

    type = ReproducingKernel{:Quadratic2D,:□,:QuinticSpline}
    integrationOrder = 4
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    cs = BenchmarkExample.sphericalCoordinate(BenchmarkExample.SphericalShell.𝑅)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
    elements["Ω"] = getCurvedElements(nodes, entities["Ω"], type, cs, integrationOrder, sp)
    elements["Γᵇ"] = getCurvedElements(nodes, entities["Γᵇ"], type, cs, integrationOrder, sp)
    elements["Γʳ"] = getCurvedElements(nodes, entities["Γʳ"], type, cs, integrationOrder, sp)
    elements["Γᵗ"] = getCurvedElements(nodes, entities["Γᵗ"], type, cs, integrationOrder, sp)
    elements["Γˡ"] = getCurvedElements(nodes, entities["Γˡ"], type, cs, integrationOrder, sp)
    elements["𝐴"] = getElements(nodes, entities["𝐴"], type, integrationOrder, sp)
    elements["𝐵"] = getElements(nodes, entities["𝐵"], type, integrationOrder, sp)

    gmsh.finalize()
    return elements, nodes
end

function import_spherical_mix(filename::String)
    gmsh.initialize()
    gmsh.open(filename)
    type = ReproducingKernel{:Quadratic2D,:□,:QuinticSpline}
    entities = getPhysicalGroups()
    nodes = get𝑿ᵢ()
    x = nodes.x
    y = nodes.y
    z = nodes.z
    sp = RegularGrid(x,y,z,n = 3,γ = 5)
    cs = BenchmarkExample.sphericalCoordinate(BenchmarkExample.SphericalShell.𝑅)
    elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()

    integrationOrder= 2
    elements["Ω"] = getCurvedPiecewiseElements(entities["Ω"], PiecewisePolynomial{:Linear2D}, cs, integrationOrder)

    integrationScheme = trilobatto3
    elements["Ωₚ"] = getCurvedPiecewiseElements(entities["Ω"], PiecewisePolynomial{:Linear2D}, cs, integrationScheme)
    elements["Ωₘ"] = getCurvedElements(nodes, entities["Ω"], type, cs, integrationScheme, sp)

    integrationScheme = lobatto3
    elements["Γₚ"] = getCurvedPiecewiseElements(entities["Γ"],PiecewisePolynomial{:Linear2D}, cs, integrationScheme,3)
    elements["Γₘ"] = getCurvedElements(nodes, entities["Γ"], type, cs, integrationScheme, sp)
    elements["Γᵇ"] = getCurvedElements(nodes, entities["Γᵇ"], type, cs, integrationScheme, sp)
    elements["Γʳ"] = getCurvedElements(nodes, entities["Γʳ"], type, cs, integrationScheme, sp)
    elements["Γᵗ"] = getCurvedElements(nodes, entities["Γᵗ"], type, cs, integrationScheme, sp)
    elements["Γˡ"] = getCurvedElements(nodes, entities["Γˡ"], type, cs, integrationScheme, sp)
    elements["Γᵇₚ"] = getElements(entities["Γᵇ"], entities["Γ"], elements["Γₚ"])
    elements["Γʳₚ"] = getElements(entities["Γʳ"], entities["Γ"], elements["Γₚ"])
    elements["Γᵗₚ"] = getElements(entities["Γᵗ"], entities["Γ"], elements["Γₚ"])
    elements["Γˡₚ"] = getElements(entities["Γˡ"], entities["Γ"], elements["Γₚ"])

    elements["𝐴"] = getElements(nodes, entities["𝐴"], type, integrationOrder, sp)
    elements["𝐵"] = getElements(nodes, entities["𝐵"], type, integrationOrder, sp)

    gmsh.finalize()
    return elements, nodes
end

prescribeForGauss = quote
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    ∂²𝗠∂x² = (0,zeros(nₘ))
    ∂²𝗠∂x∂y = (0,zeros(nₘ))
    ∂²𝗠∂y² = (0,zeros(nₘ))
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
    push!(elements["Ω"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠, :∂²𝝭∂x²=>:𝑠, :∂²𝝭∂x∂y=>:𝑠, :∂²𝝭∂y²=>:𝑠)
    push!(elements["Ω"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂²𝗠∂x²=>∂²𝗠∂x², :∂²𝗠∂y²=>∂²𝗠∂y², :∂²𝗠∂x∂y=>∂²𝗠∂x∂y)
end
prescribeForMix = quote
    𝗠 = (0,zeros(nₘ))
    ∂𝗠∂x = (0,zeros(nₘ))
    ∂𝗠∂y = (0,zeros(nₘ))
    ∂²𝗠∂x² = (0,zeros(nₘ))
    ∂²𝗠∂x∂y = (0,zeros(nₘ))
    ∂²𝗠∂y² = (0,zeros(nₘ))
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

    prescribe!(elements["Γₚ"],:𝒂₃₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γₚ"],:𝒂₃₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γₚ"],:𝒂₃₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γₚ"],:Γ¹₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ¹₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ¹₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ²₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ²₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γₚ"],:Γ²₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))

    push!(elements["Ωₘ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Ωₘ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γₘ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γₘ"], :𝗠=>𝗠, :∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)

    push!(elements["Ω"], :𝝭=>:𝑠)
    push!(elements["Ωₚ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠, :∂²𝝭∂x²=>:𝑠, :∂²𝝭∂x∂y=>:𝑠, :∂²𝝭∂y²=>:𝑠)
    push!(elements["Γₚ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)

    push!(elements["Γˡ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γʳ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γˡ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γʳ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
end

prescribeForPenalty = quote
    prescribe!(elements["Γˡ"],:𝒂₃₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γˡ"],:𝒂₃₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γˡ"],:𝒂₃₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γˡ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₁₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₂₂=>(ξ¹,ξ²,ξ³)->1.0)
    prescribe!(elements["Γˡ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₃₃=>(ξ¹,ξ²,ξ³)->0.0)

    prescribe!(elements["Γʳ"],:𝒂₃₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γʳ"],:𝒂₃₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γʳ"],:𝒂₃₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γʳ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:n₁₁=>(ξ¹,ξ²,ξ³)->1.0)
    prescribe!(elements["Γʳ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:n₂₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:n₃₃=>(ξ¹,ξ²,ξ³)->0.0)

    push!(elements["Γˡ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γʳ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠)
    push!(elements["Γˡ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
    push!(elements["Γʳ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y)
end
prescribeForNitsche = quote
    ∂³𝗠∂x³ = (0,zeros(nₘ))
    ∂³𝗠∂x²∂y = (0,zeros(nₘ))
    ∂³𝗠∂x∂y² = (0,zeros(nₘ))
    ∂³𝗠∂y³ = (0,zeros(nₘ))
    prescribe!(elements["Γʳ"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:a²²=>(ξ¹,ξ²,ξ³)->cs.a²²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:a¹²=>(ξ¹,ξ²,ξ³)->cs.a¹²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:b₁₁=>(ξ¹,ξ²,ξ³)->cs.b₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:b₂₂=>(ξ¹,ξ²,ξ³)->cs.b₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:b₁₂=>(ξ¹,ξ²,ξ³)->cs.b₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:𝒂₁₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γʳ"],:𝒂₁₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γʳ"],:𝒂₁₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γʳ"],:𝒂₂₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γʳ"],:𝒂₂₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γʳ"],:𝒂₂₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γʳ"],:𝒂₃₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γʳ"],:𝒂₃₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γʳ"],:𝒂₃₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γʳ"],:𝒂¹₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γʳ"],:𝒂¹₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γʳ"],:𝒂¹₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂¹(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γʳ"],:𝒂²₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γʳ"],:𝒂²₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γʳ"],:𝒂²₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂²(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γʳ"],:𝒂³₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂³(Vec{3}((ξ¹,ξ²,ξ³)))[1])
    prescribe!(elements["Γʳ"],:𝒂³₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂³(Vec{3}((ξ¹,ξ²,ξ³)))[2])
    prescribe!(elements["Γʳ"],:𝒂³₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂³(Vec{3}((ξ¹,ξ²,ξ³)))[3])
    prescribe!(elements["Γʳ"],:Γ¹₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ¹₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ¹₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ²₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ²₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ²₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ¹₁₁₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ¹₁₁₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ¹₁₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ¹₁₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ¹₂₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ¹₂₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ¹₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ²₁₁₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ²₁₁₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ²₁₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ²₁₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ²₂₂₁=>(ξ¹,ξ²,ξ³)->cs.∂₁Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γʳ"],:Γ²₂₂₂=>(ξ¹,ξ²,ξ³)->cs.∂₂Γ²₂₂(Vec{3}((ξ¹,ξ²,ξ³))))

    prescribe!(elements["Γˡ"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:a²²=>(ξ¹,ξ²,ξ³)->cs.a²²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:a¹²=>(ξ¹,ξ²,ξ³)->cs.a¹²(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:b₁₁=>(ξ¹,ξ²,ξ³)->cs.b₁₁(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:b₂₂=>(ξ¹,ξ²,ξ³)->cs.b₂₂(Vec{3}((ξ¹,ξ²,ξ³))))
    prescribe!(elements["Γˡ"],:b₁₂=>(ξ¹,ξ²,ξ³)->cs.b₁₂(Vec{3}((ξ¹,ξ²,ξ³))))
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

    push!(elements["Γʳ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠, :∂²𝝭∂x²=>:𝑠, :∂²𝝭∂x∂y=>:𝑠, :∂²𝝭∂y²=>:𝑠, :∂³𝝭∂x³=>:𝑠, :∂³𝝭∂x²∂y=>:𝑠, :∂³𝝭∂x∂y²=>:𝑠, :∂³𝝭∂y³=>:𝑠)
    push!(elements["Γˡ"], :𝝭=>:𝑠, :∂𝝭∂x=>:𝑠, :∂𝝭∂y=>:𝑠, :∂²𝝭∂x²=>:𝑠, :∂²𝝭∂x∂y=>:𝑠, :∂²𝝭∂y²=>:𝑠, :∂³𝝭∂x³=>:𝑠, :∂³𝝭∂x²∂y=>:𝑠, :∂³𝝭∂x∂y²=>:𝑠, :∂³𝝭∂y³=>:𝑠)
    push!(elements["Γʳ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂²𝗠∂x²=>∂²𝗠∂x², :∂²𝗠∂y²=>∂²𝗠∂y², :∂²𝗠∂x∂y=>∂²𝗠∂x∂y, :∂³𝗠∂x³=>∂³𝗠∂x³, :∂³𝗠∂x²∂y=>∂³𝗠∂x²∂y, :∂³𝗠∂x∂y²=>∂³𝗠∂x∂y², :∂³𝗠∂y³=>∂³𝗠∂y³)
    push!(elements["Γˡ"], :𝗠=>𝗠,:∂𝗠∂x=>∂𝗠∂x, :∂𝗠∂y=>∂𝗠∂y, :∂²𝗠∂x²=>∂²𝗠∂x², :∂²𝗠∂y²=>∂²𝗠∂y², :∂²𝗠∂x∂y=>∂²𝗠∂x∂y, :∂³𝗠∂x³=>∂³𝗠∂x³, :∂³𝗠∂x²∂y=>∂³𝗠∂x²∂y, :∂³𝗠∂x∂y²=>∂³𝗠∂x∂y², :∂³𝗠∂y³=>∂³𝗠∂y³)
end

prescribeVariables = quote
    prescribe!(elements["Γʳ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:n₁₁=>(ξ¹,ξ²,ξ³)->1.0)
    prescribe!(elements["Γʳ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:n₂₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γʳ"],:n₃₃=>(ξ¹,ξ²,ξ³)->0.0)

    prescribe!(elements["Γˡ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₁₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₂₂=>(ξ¹,ξ²,ξ³)->1.0)
    prescribe!(elements["Γˡ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"],:n₃₃=>(ξ¹,ξ²,ξ³)->0.0)

    prescribe!(elements["Γʳ"], :θ=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["Γˡ"], :θ=>(ξ¹,ξ²,ξ³)->0.0)

    push!(elements["Γʳ"], :∂ₙu₁=>:𝐺, :∂ₙu₂=>:𝐺, :∂ₙu₃=>:𝐺)
    push!(elements["Γˡ"], :∂ₙu₁=>:𝐺, :∂ₙu₂=>:𝐺, :∂ₙu₃=>:𝐺)
    for a in elements["Γʳ"]∪elements["Γˡ"]
        for ξ in a.𝓖
            n¹ = ξ.n¹
            n² = ξ.n²
            ξ.∂ₙu₁ = 0.0
            ξ.∂ₙu₂ = 0.0
            ξ.∂ₙu₃ = 0.0
        end
    end
    prescribe!(elements["𝐴"],:t₁=>(ξ¹,ξ²,ξ³)->2.0)
    prescribe!(elements["𝐴"],:t₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["𝐴"],:t₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["𝐴"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["𝐴"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["𝐴"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["𝐴"],:n₁₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["𝐴"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["𝐴"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["𝐴"],:n₂₂=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["𝐴"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["𝐴"],:n₃₃=>(ξ¹,ξ²,ξ³)->1.0)
    prescribe!(elements["𝐵"],:t₁=>(ξ¹,ξ²,ξ³)->0.0)
    prescribe!(elements["𝐵"],:t₂=>(ξ¹,ξ²,ξ³)->-2.0)
    prescribe!(elements["𝐵"],:t₃=>(ξ¹,ξ²,ξ³)->0.0)

    push!(elements["𝐴"], :𝝭=>:𝑠)
    push!(elements["𝐴"], :𝗠=>𝗠)
    push!(elements["𝐵"], :𝝭=>:𝑠)
    push!(elements["𝐵"], :𝗠=>𝗠)
end
