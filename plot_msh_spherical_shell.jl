
using ApproxOperator, GLMakie

import BenchmarkExample: BenchmarkExample
import Gmsh: gmsh

𝑅 = BenchmarkExample.SphericalShell.𝑅

ndiv = 40
gmsh.initialize()
gmsh.open("./msh/sphericalshell_"*string(ndiv)*".msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()
ξ¹ = nodes.x
ξ² = nodes.y
elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes,entities["Ω"])
elements["Γ¹"] = getElements(nodes,entities["Γᵇ"])
elements["Γ²"] = getElements(nodes,entities["Γʳ"])
elements["Γ³"] = getElements(nodes,entities["Γᵗ"])
elements["Γ⁴"] = getElements(nodes,entities["Γˡ"])
elements["∂Ω"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]

gmsh.finalize()

f = Figure()

# axis
ax = Axis3(f[1, 1], perspectiveness = 0.8, aspect = :data, azimuth = 0.25*pi, elevation = 0.2*pi, xlabel = " ", ylabel = " ", zlabel = " ", xticksvisible = false,xticklabelsvisible=false, yticksvisible = false, yticklabelsvisible=false, zticksvisible = false, zticklabelsvisible=false, protrusions = (0.,0.,0.,0.))
hidespines!(ax)
# hidedecorations!(ax)

# x = ξ¹
# y = ξ²
# z = zeros(length(ξ¹))
# ps = Point3f.(x,y,z)
# scatter!(ps, 
#     marker=:circle,
#     markersize = 5,
#     color = :black
# )

# # elements
# for elm in elements["Ω"]
#     x = [x.x for x in elm.𝓒[[1,2,3,1]]]
#     y = [x.y for x in elm.𝓒[[1,2,3,1]]]
#     lines!(x,y,linestyle = :dash, linewidth = 0.5, color = :black)
# end

# # boundaries
# for elm in elements["∂Ω"]
#     x = [x.x for x in elm.𝓒]
#     y = [x.y for x in elm.𝓒]
#     lines!(x,y,linewidth = 1.5, color = :black)
# end

# nodes
x = [𝑅*cos(ξ/𝑅)*cos(η/𝑅) for (ξ,η) in zip(ξ¹,ξ²)]
y = [𝑅*sin(ξ/𝑅)*cos(η/𝑅) for (ξ,η) in zip(ξ¹,ξ²)]
z = [𝑅*sin(ξ/𝑅) for ξ in ξ²]
ps = Point3f.(x,y,z)
scatter!(ps, 
    marker=:circle,
    markersize = 3,
    color = :black
)

# elements
for elm in elements["Ω"]
    x = [𝑅*cos(x.x/𝑅)*cos(x.y/𝑅) for x in elm.𝓒[[1,2,3,1]]]
    y = [𝑅*sin(x.x/𝑅)*cos(x.y/𝑅) for x in elm.𝓒[[1,2,3,1]]]
    z = [𝑅*sin(x.y/𝑅) for x in elm.𝓒[[1,2,3,1]]]
    lines!(x,y,z,linestyle = :dash, linewidth = 0.3, color = :black)
end

# boundaries
for elm in elements["∂Ω"]
    ξ¹ = [x.x for x in elm.𝓒]
    ξ² = [x.y for x in elm.𝓒]
    x = [𝑅*cos(ξ/𝑅)*cos(η/𝑅) for (ξ,η) in zip(ξ¹,ξ²)]
    y = [𝑅*sin(ξ/𝑅)*cos(η/𝑅) for (ξ,η) in zip(ξ¹,ξ²)]
    z = [𝑅*sin(ξ/𝑅) for ξ in ξ²]
    lines!(x,y,z,linewidth = 1.5, color = :black)
end

save("./png/sphericalshell_"*string(ndiv)*".png",f, px_per_unit = 10.0)
f