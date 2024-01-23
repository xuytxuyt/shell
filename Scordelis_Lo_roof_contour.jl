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
# gmsh.finalize()

type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
𝗠 = zeros(21)

d₁, d₂, d₃ = load("jld/scordelislo_gauss_"*string(ndiv)*".jld")
# d₁, d₂, d₃ = load("jld/scordelislo_mix_"*string(ndiv)*".jld")
push!(nodes,:d₁=>d₁[2],:d₂=>d₂[2],:d₃=>d₃[2])

ind = 11
xs = zeros(ind)
x₂ = zeros(ind)
ys = zeros(ind)
y₂ = zeros(ind)
zs = zeros(ind,ind)
z₂ = zeros(ind,ind)
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
        xs[I] = 𝑅*sin(ξ¹/𝑅) 
        x₂[I] = 𝑅*sin(ξ¹/𝑅-𝜃)
        ys[J] = ξ² + α*u₂
        y₂[J] = 50-ξ² - α*u₂
        zs[I,J] = 𝑅*cos(ξ¹/𝑅) + α*u₃
        z₂[I,J] = 𝑅*cos(𝜃-ξ¹/𝑅) + α*u₃
        color[I,J] =u₃
    end
end

xs1 = zeros(ind)
ys1 = zeros(ind)
zs1 = zeros(ind,ind)
color1 = zeros(ind)
for i in 1:ind
    xs1[i] = -xs[ind-i+1]
    for j in 1:ind
        ys1[j] = ys[j]
        zs1[i,j] = zs[ind-i+1,j]
        color1[i,j] = color1[ind-i+1,j]
    end
end

fig = Figure()
ax = Axis3(fig[1, 1], aspect = (1, 1, 0.1),azimuth = 0.3*pi, elevation = 0.2*pi)

# hidespines!(ax)
# hidedecorations!(ax)
s = surface!(ax,xs,ys,zs, color=color, colormap=:redsblues)
surface!(ax,x₂,ys,z₂, color=color, colormap=:redsblues)
surface!(ax,xs,y₂,zs, color=color, colormap=:redsblues)
surface!(ax,x₂,y₂,zs, color=color, colormap=:redsblues)
lines!([Point(16.0696902421634,0,19.15111107412126),Point(16.0696902421634,50,19.15111107412126)],color=:black)
lines!([Point(-16.0696902421634,0,19.15111107412126),Point(-16.0696902421634,50,19.15111107412126)],color=:black)
Colorbar(fig[2, 1], s, vertical = false)

fig