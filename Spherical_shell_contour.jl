using ApproxOperator, JLD, GLMakie

import Gmsh: gmsh
import BenchmarkExample: BenchmarkExample

𝑅 = BenchmarkExample.SphericalShell.𝑅
F = BenchmarkExample.SphericalShell.𝐹
E = BenchmarkExample.SphericalShell.𝐸
ν = BenchmarkExample.SphericalShell.𝜈
h = BenchmarkExample.SphericalShell.ℎ
𝜃₁ =  BenchmarkExample.SphericalShell.𝜃₁
𝜃₂ =  BenchmarkExample.SphericalShell.𝜃₂
cs = BenchmarkExample.sphericalCoordinate(𝑅)

ndiv = 32
α = 1e1

## import nodes
gmsh.initialize()
gmsh.open("msh/sphericalshell_"*string(ndiv)*".msh")
nodes = get𝑿ᵢ()
x = nodes.x
y = nodes.y
z = nodes.z
sp = RegularGrid(x,y,z,n = 3,γ = 5)
nₚ = length(nodes)
s = 2.5*𝑅*𝜃/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
# gmsh.finalize()

type = ReproducingKernel{:Quadratic2D,:□,:QuinticSpline}
𝗠 = zeros(21)
# type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
# 𝗠 = zeros(55)

ds = Dict(load("jld/spherical_shell_mix_hr_"*string(ndiv)*".jld"))
push!(nodes,:d₁=>ds["d₁"],:d₂=>ds["d₂"],:d₃=>ds["d₃"])

ind = 21
xs1 = zeros(ind,ind)
ys1 = zeros(ind,ind)
zs1 = zeros(ind,ind)
cs1 = zeros(ind,ind)
xs2 = zeros(ind,ind)
ys2 = zeros(ind,ind)
zs2 = zeros(ind,ind)
cs2 = zeros(ind,ind)
xs3 = zeros(ind,ind)
ys3 = zeros(ind,ind)
zs3 = zeros(ind,ind)
cs3 = zeros(ind,ind)
xs4 = zeros(ind,ind)
ys4 = zeros(ind,ind)
zs4 = zeros(ind,ind)
cs4 = zeros(ind,ind)
for (I,ξ¹) in enumerate(LinRange(0.0, 𝑅*𝜃₂, ind))
    for (J,ξ²) in enumerate(LinRange(0.0, 𝑅*𝜃₁, ind))
        indices = sp(ξ¹,ξ²,0.0)
        N = zeros(length(indices))
        data = Dict([:x=>(2,[ξ¹]),:y=>(2,[ξ²]),:z=>(2,[0.0]),:𝝭=>(4,N),:𝗠=>(0,𝗠)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set𝝭!(ap)
        u₁ = 0.0
        u₂ = 0.0
        u₃ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            u₃ += N[i]*xᵢ.d₃
        end
        x_ = Vec{3}((ξ¹,ξ²,0.0))
        𝒂₃ = cs.𝒂₃(x_)
        u = Vec{3}((u₁,u₂,u₃)) ⋅ 𝒂₃
        xs1[I,J] = 𝑅*cos(ξ¹/𝑅)*cos(ξ²/𝑅) + α*u₁
        ys1[I,J] = 𝑅*sin(ξ¹/𝑅)*cos(ξ²/𝑅) + α*u₂
        zs1[I,J] = 𝑅*sin(ξ²/𝑅) + α*u₃
        cs1[I,J] = u
    end
end

for I in 1:ind
    for J in 1:ind
        xs2[I,J] = -xs1[ind-I+1,J]
        ys2[I,J] = ys1[ind-I+1,J]
        zs2[I,J] = zs1[ind-I+1,J]
        cs2[I,J] = cs1[ind-I+1,J]
        xs3[I,J] = -xs1[ind-I+1,ind-J+1]
        ys3[I,J] = 𝐿-ys1[ind-I+1,ind-J+1]
        zs3[I,J] = zs1[ind-I+1,ind-J+1]
        cs3[I,J] = cs1[ind-I+1,ind-J+1]
        xs4[I,J] = xs1[I,ind-J+1]
        ys4[I,J] = 𝐿-ys1[I,ind-J+1]
        zs4[I,J] = zs1[I,ind-J+1]
        cs4[I,J] = cs1[I,ind-J+1]
    end
end

fig = Figure()
ax = Axis3(fig[1, 1], perspectiveness = 0.8, aspect = :equal, azimuth = 0.25*pi, elevation = 0.2*pi, xlabel = " ", ylabel = " ", zlabel = " ", xticksvisible = false,xticklabelsvisible=false, yticksvisible = false, yticklabelsvisible=false, zticksvisible = false, zticklabelsvisible=false, protrusions = (0.,0.,0.,0.))
# ax = Axis3(fig[1, 1], perspectiveness = 1, aspect = :data, azimuth = -0.4*pi, elevation = 0.2*pi)

# hidespines!(ax)
# hidedecorations!(ax)
s = surface!(ax,xs1,ys1,zs1, color=cs1, colormap=:lightrainbow)
# surface!(ax,xs2,ys2,zs2, color=cs2, colormap=:redsblues)
# surface!(ax,xs3,ys3,zs3, color=cs3, colormap=:redsblues)
# surface!(ax,xs4,ys4,zs4, color=cs4, colormap=:redsblues)
# lines!(ax,xl₁,yl₁,zl₁,color=:black)
# lines!(ax,xl₂,yl₂,zl₂,color=:black)
# lines!(ax,xl₃,yl₃,zl₃,color=:black)
# lines!(ax,xl₄,yl₄,zl₄,color=:black)
# lines!(ax,xl₅,yl₅,zl₅,color=:black)
# lines!(ax,xl₆,yl₆,zl₆,color=:black)
# lines!(ax,xl₇,yl₇,zl₇,color=:black)
# lines!(ax,xl₈,yl₈,zl₈,color=:black)
Colorbar(fig[2, 1], s, vertical = false)

fig
