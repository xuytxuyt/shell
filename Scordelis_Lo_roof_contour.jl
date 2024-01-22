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

ndiv = 11
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
gmsh.finalize()

type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
𝗠 = zeros(21)

# d₁, d₂, d₃ = load("jld/scordelislo_gauss"*string(ndiv)*".jld")
d₁, d₂, d₃ = load("jld/scordelislo_mix_"*string(ndiv)*".jld")
push!(nodes,:d₁=>d₁[2],:d₂=>d₂[2],:d₃=>d₃[2])

ind = 10
xs = zeros(ind)
ys = zeros(ind)
zs = zeros(ind,ind)
color = zeros(ind,ind)
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
        xs[I] = 𝑅*sin(ξ¹/𝑅) + α*u₁
        ys[J] = ξ² + α*u₂
        zs[I,J] = 𝑅*cos(ξ¹/𝑅) + α*u₃
        color[I,J] = u₃
    end
end

fig = Figure()

ax = Axis3(fig[1, 1])

hidespines!(ax)
hidedecorations!(ax)
# lines!([Point(0, 0, 25), Point(17.45329251994329,0,25), Point(17.45329251994329, 25, 25), Point(0, 25, 25), Point(0, 0, 25)],color=:black)
s = surface!(ax,xs,ys,zs, color=color, colormap=:redsblues)
Colorbar(fig[2, 1], s, vertical = false)
fig