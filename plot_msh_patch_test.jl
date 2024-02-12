
using ApproxOperator, GLMakie

import BenchmarkExample: BenchmarkExample
import Gmsh: gmsh

𝑅 = BenchmarkExample.PatchTestThinShell.𝑅

filename = "patchtest"
gmsh.initialize()
gmsh.open("./msh/"*filename*".msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()
ξ¹ = nodes.x
ξ² = nodes.y
elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes,entities["Ω"])
elements["Γ¹"] = getElements(nodes,entities["Γ¹"])
elements["Γ²"] = getElements(nodes,entities["Γ²"])
elements["Γ³"] = getElements(nodes,entities["Γ³"])
elements["Γ⁴"] = getElements(nodes,entities["Γ⁴"])
elements["∂Ω"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]

gmsh.finalize()

f = Figure()

# axis
ax = Axis3(f[1, 1], perspectiveness = 0.8, aspect = :equal, azimuth = -0.25*pi, elevation = 0.2*pi, xlabel = " ", ylabel = " ", zlabel = " ", xticksvisible = false,xticklabelsvisible=false, yticksvisible = false, yticklabelsvisible=false, zticksvisible = false, zticklabelsvisible=false, protrusions = (0.,0.,0.,0.))
hidespines!(ax)
# hidedecorations!(ax)

x = ξ¹
y = ξ²
z = zeros(length(ξ¹))
ps = Point3f.(x,y,z)
scatter!(ps, 
    marker=:circle,
    markersize = 5,
    color = :black
)

# elements
for elm in elements["Ω"]
    x = [x.x for x in elm.𝓒[[1,2,3,1]]]
    y = [x.y for x in elm.𝓒[[1,2,3,1]]]
    lines!(x,y,linestyle = :dash, linewidth = 0.5, color = :black)
end

# boundaries
for elm in elements["∂Ω"]
    x = [x.x for x in elm.𝓒]
    y = [x.y for x in elm.𝓒]
    lines!(x,y,linewidth = 1.5, color = :black)
end

# # nodes
# x = [𝑅*sin(ξ/𝑅) for ξ in ξ¹]
# y = ξ²
# z = [𝑅*cos(ξ/𝑅) for ξ in ξ¹]
# ps = Point3f.(x,y,z)
# scatter!(ps, 
#     marker=:circle,
#     markersize = 5,
#     color = :black
# )

# # elements
# for elm in elements["Ω"]
#     x = [𝑅*sin(x.x/𝑅) for x in elm.𝓒[[1,2,3,1]]]
#     y = [x.y for x in elm.𝓒[[1,2,3,1]]]
#     z = [𝑅*cos(x.x/𝑅) for x in elm.𝓒[[1,2,3,1]]]
#     lines!(x,y,z,linestyle = :dash, linewidth = 0.5, color = :black)
# end

# # boundaries
# for elm in elements["∂Ω"]
#     ξ¹ = [x.x for x in elm.𝓒]
#     ξ² = [x.y for x in elm.𝓒]
#     x = [𝑅*sin(ξ/𝑅) for ξ in ξ¹]
#     y = ξ²
#     z = [𝑅*cos(ξ/𝑅) for ξ in ξ¹]
#     lines!(x,y,z,linewidth = 1.5, color = :black)
# end

# arc(Point3f(0), 0.1,0.0, 1.0)

# arrow
# arrows!([0.,0.,0.],[0.,0.,0.],[0.,0.,0.],[1,0.,0.],[0.,1,0.],[0.,0.,1], arrowsize = 0.03, lengthscale = 0.5, linewidth = 0.005)
# lines!([0.,𝑅*sin(1.0/𝑅)],[0.,0.],[0.,𝑅*cos(1.0/𝑅)],linewidth = 0.5, color = :black)
# lines!([0.,0.],[0.,0.],[0.,𝑅*cos(0.0/𝑅)],linewidth = 0.5, color = :black)
# save("./png/"*filename*"_curved_msh.png",f, px_per_unit = 10.0)
save("./png/"*filename*"_flat_msh.png",f, px_per_unit = 10.0)
f