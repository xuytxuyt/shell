using ApproxOperator, JLD, GLMakie

import Gmsh: gmsh
import BenchmarkExample: BenchmarkExample

𝑅 = BenchmarkExample.ScordelisLoRoof.𝑅
𝐿 = BenchmarkExample.ScordelisLoRoof.𝐿
𝜃 = BenchmarkExample.ScordelisLoRoof.𝜃
b₃ = BenchmarkExample.ScordelisLoRoof.𝑞
E = BenchmarkExample.ScordelisLoRoof.𝐸
ν = BenchmarkExample.ScordelisLoRoof.𝜈
h = BenchmarkExample.ScordelisLoRoof.ℎ
cs = BenchmarkExample.cylindricalCoordinate(𝑅)

ndiv = 16
α = 10.0

## import nodes
gmsh.initialize()
gmsh.open("msh/scordelislo_"*string(ndiv)*".msh")
nodes = get𝑿ᵢ()
x = nodes.x
y = nodes.y
z = nodes.z
sp = RegularGrid(x,y,z,n = 3,γ = 5)
nₚ = length(nodes)
s = 2.5*𝐿/2/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
# gmsh.finalize()

type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
𝗠 = zeros(21)

d₁, d₂, d₃ = load("jld/scordelislo_gauss_"*string(ndiv)*".jld")
# d₁, d₂, d₃ = load("jld/scordelislo_mix_"*string(ndiv)*".jld")
push!(nodes,:d₁=>d₁[2],:d₂=>d₂[2],:d₃=>d₃[2])

ind = 11
xs1 = zeros(ind)
ys1 = zeros(ind)
zs1 = zeros(ind,ind)
cs1 = zeros(ind,ind)
xs2 = zeros(ind)
ys2 = zeros(ind)
zs2 = zeros(ind,ind)
cs2 = zeros(ind,ind)
xs3 = zeros(ind)
ys3 = zeros(ind)
zs3 = zeros(ind,ind)
cs3 = zeros(ind,ind)
xs4 = zeros(ind)
ys4 = zeros(ind)
zs4 = zeros(ind,ind)
cs4 = zeros(ind,ind)
xl₁ = zeros(ind)
xl₂ = zeros(ind)
xl₃ = zeros(ind)
xl₄ = zeros(ind)
xl₅ = zeros(ind)
xl₆ = zeros(ind)
xl₇ = zeros(ind)
xl₈ = zeros(ind)
yl₁ = zeros(ind)
yl₂ = zeros(ind)
yl₃ = zeros(ind)
yl₄ = zeros(ind)
yl₅ = zeros(ind)
yl₆ = zeros(ind)
yl₇ = zeros(ind)
yl₈ = zeros(ind)
zl₁ = zeros(ind)
zl₂ = zeros(ind)
zl₃ = zeros(ind)
zl₄ = zeros(ind)
zl₅ = zeros(ind)
zl₆ = zeros(ind)
zl₇ = zeros(ind)
zl₈ = zeros(ind)
for (I,ξ¹) in enumerate(LinRange(0.0, 𝜃*𝑅, ind))
    for (J,ξ²) in enumerate(LinRange(0.0, 𝐿/2, ind))
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
        xs1[I] = 𝑅*sin(ξ¹/𝑅)
        ys1[J] = ξ²
        zs1[I,J] = 𝑅*cos(ξ¹/𝑅)
        cs1[I,J] = u₂
        xl₂[J] = 𝑅*sin(𝜃)
        yl₂[J] = ξ²
        zl₂[J] = 𝑅*cos(𝜃)
        xl₄[J] =-𝑅*sin(𝜃)
        yl₄[J] = ξ²
        zl₄[J] = 𝑅*cos(𝜃)
        xl₆[J] = 𝑅*sin(𝜃)
        yl₆[J] = ξ²+𝐿/2
        zl₆[J] = 𝑅*cos(𝜃)
        xl₈[J] =-𝑅*sin(𝜃)
        yl₈[J] = ξ²+𝐿/2
        zl₈[J] = 𝑅*cos(𝜃)
    end
    xl₁[I] = 𝑅*sin(ξ¹/𝑅)
    yl₁[I] = 0.0
    zl₁[I] = 𝑅*cos(ξ¹/𝑅)
    xl₃[I] = 𝑅*sin(ξ¹/𝑅)
    yl₃[I] = 𝐿
    zl₃[I] = 𝑅*cos(ξ¹/𝑅)
    xl₅[I] =-𝑅*sin(ξ¹/𝑅)
    yl₅[I] = 0.0
    zl₅[I] = 𝑅*cos(ξ¹/𝑅)
    xl₇[I] =-𝑅*sin(ξ¹/𝑅)
    yl₇[I] = 𝐿
    zl₇[I] = 𝑅*cos(ξ¹/𝑅)
end

for I in 1:ind
    for J in 1:ind
        xs2[I] = -xs1[ind-I+1]
        ys2[J] = ys1[J]
        zs2[I,J] = zs1[ind-I+1,J]
        cs2[I,J] = cs1[ind-I+1,J]
        xs3[I] = -xs1[ind-I+1]
        ys3[J] = 𝐿-ys1[ind-J+1]
        zs3[I,J] = zs1[ind-I+1,ind-J+1]
        cs3[I,J] = cs1[ind-I+1,ind-J+1]
        xs4[I] = xs1[I]
        ys4[J] = 𝐿-ys1[ind-J+1]
        zs4[I,J] = zs1[I,ind-J+1]
        cs4[I,J] = cs1[I,ind-J+1]
    end
end
# xs1 = zeros(ind)
# ys1 = zeros(ind)
# zs1 = zeros(ind,ind)
# color1 = zeros(ind)
# for i in 1:ind
#     xs1[i] = -xs[ind-i+1]
#     for j in 1:ind
#         ys1[j] = ys[j]
#         zs1[i,j] = zs[ind-i+1,j]
#         color1[i,j] = color1[ind-i+1,j]
#     end
# end

fig = Figure()
ax = Axis3(fig[1, 1], perspectiveness = 1, aspect = :data, azimuth = -0.4*pi, elevation = 0.2*pi)

hidespines!(ax)
hidedecorations!(ax)
s = surface!(ax,xs1,ys1,zs1, color=cs1, colormap=:redsblues)
surface!(ax,xs2,ys2,zs2, color=cs2, colormap=:redsblues)
surface!(ax,xs3,ys3,zs3, color=cs3, colormap=:redsblues)
surface!(ax,xs4,ys4,zs4, color=cs4, colormap=:redsblues)
lines!(ax,xl₁,yl₁,zl₁,color=:black)
lines!(ax,xl₂,yl₂,zl₂,color=:black)
lines!(ax,xl₃,yl₃,zl₃,color=:black)
lines!(ax,xl₄,yl₄,zl₄,color=:black)
lines!(ax,xl₅,yl₅,zl₅,color=:black)
lines!(ax,xl₆,yl₆,zl₆,color=:black)
lines!(ax,xl₇,yl₇,zl₇,color=:black)
lines!(ax,xl₈,yl₈,zl₈,color=:black)
Colorbar(fig[2, 1], s, vertical = false)

fig
