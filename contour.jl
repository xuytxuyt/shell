using ApproxOperator, JLD, GLMakie, Pardiso

import Gmsh: gmsh
import BenchmarkExample: BenchmarkExample

# 𝑅 = BenchmarkExample.ScordelisLoRoof.𝑅
# 𝐿 = BenchmarkExample.ScordelisLoRoof.𝐿
# 𝜃 = BenchmarkExample.ScordelisLoRoof.𝜃
# b₃ = BenchmarkExample.ScordelisLoRoof.𝑞
# E = BenchmarkExample.ScordelisLoRoof.𝐸
# ν = BenchmarkExample.ScordelisLoRoof.𝜈
# h = BenchmarkExample.ScordelisLoRoof.ℎ
# cs = BenchmarkExample.cylindricalCoordinate(𝑅)

ndiv = 12
α = 1e2

## import nodes
gmsh.initialize()
gmsh.open("msh/TADAS dampers.msh")
nodes = get𝑿ᵢ()
x = nodes.x
y = nodes.y
z = nodes.z
sp = RegularGrid(x,y,z,n = 3,γ = 5)
nₚ = length(nodes)
s = 3.5*120/ndiv*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
# gmsh.finalize()

# type = ReproducingKernel{:Quadratic2D,:□,:CubicSpline}
# 𝗠 = zeros(21)
type = ReproducingKernel{:Cubic2D,:□,:QuinticSpline}
𝗠 = zeros(55)

ds = Dict(load("jld/TADAS_hr_12.jld"))
push!(nodes,:d=>ds["d"])

ind = 21
xs = zeros(ind,ind)
ys = zeros(ind,ind)
zs = zeros(ind,ind)
cs = zeros(ind,ind)
xs1 = zeros(ind,ind)
ys1 = zeros(ind,ind)
zs1 = zeros(ind,ind)
cs1 = zeros(ind,ind)
for (I,ξ¹) in enumerate(LinRange(0.0,120, ind))
    for (J,ξ²) in enumerate(LinRange(0.0, 30, ind))
        indices = sp(ξ¹,ξ²,0.0)
        N = zeros(length(indices))
        data = Dict([:x=>(2,[ξ¹]),:y=>(2,[ξ²]),:z=>(2,[0.0]),:𝝭=>(4,N),:𝗠=>(0,𝗠)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set𝝭!(ap)
        u₃ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
        end
        if ξ² <= 30
           xs1[I,J] = ξ¹
           ys1[I,J] = ξ²
           zs1[I,J] = α*u₃
           cs1[I,J] = u₃
        end
    end
end
for (I,ξ¹) in enumerate(LinRange(30,90, ind))
    for (J,ξ²) in enumerate(LinRange(30, 80, ind))
        indices = sp(ξ¹,ξ²,0.0)
        N = zeros(length(indices))
        data = Dict([:x=>(2,[ξ¹]),:y=>(2,[ξ²]),:z=>(2,[0.0]),:𝝭=>(4,N),:𝗠=>(0,𝗠)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set𝝭!(ap)
        u₃ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
        end
        if ξ¹ <=45 && ξ²<=10/3*ξ¹-70
           xs[I,J] = ξ¹
           ys[I,J] = ξ²
           zs[I,J] = α*u₃
           cs[I,J] = u₃
        elseif ξ¹ >=75 && ξ²<=-10/3*ξ¹+330
            xs[I,J] = ξ¹
            ys[I,J] = ξ²
            zs[I,J] = α*u₃
            cs[I,J] = u₃
        elseif ξ¹ >=45 && ξ¹ <=75
            xs[I,J] = ξ¹
            ys[I,J] = ξ²
            zs[I,J] = α*u₃
            cs[I,J] = u₃
        elseif ξ¹ <=45
            xs[I,J] = ξ¹
            ys[I,J] = 10/3*ξ¹-70
            zs[I,J] = α*u₃
            cs[I,J] = u₃
        elseif ξ¹ >=75
            xs[I,J] = ξ¹
            ys[I,J] = -10/3*ξ¹+330
            zs[I,J] = α*u₃
            cs[I,J] = u₃
        end
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

# hidespines!(ax)
# hidedecorations!(ax)
s = surface!(ax,xs,ys,zs, color=cs, colormap=:redsblues)
s = surface!(ax,xs1,ys1,zs1, color=cs1, colormap=:redsblues)

Colorbar(fig[2, 1], s, vertical = false)

fig
