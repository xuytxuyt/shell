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
α = 2e2

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
type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
𝗠 = zeros(55)
h = 5
ds = Dict(load("jld/TADAS_hr_12.jld"))
push!(nodes,:d=>ds["d"])

ind = 21
xs1 = zeros(ind,ind)
ys1 = zeros(ind,ind)
zs1 = zeros(ind,ind)
cs1 = zeros(ind,ind)
ys1s = zeros(ind,ind)
zs1s = zeros(ind,ind)
ys1x = zeros(ind,ind)
zs1x = zeros(ind,ind)
xs2 = zeros(ind,ind)
ys2 = zeros(ind,ind)
zs2 = zeros(ind,ind)
cs2 = zeros(ind,ind)
ys2x = zeros(ind,ind)
zs2x = zeros(ind,ind)
ys2s = zeros(ind,ind)
zs2s = zeros(ind,ind)
xs3 = zeros(ind,ind)
ys3 = zeros(ind,ind)
zs3 = zeros(ind,ind)
cs3 = zeros(ind,ind)
ys3s = zeros(ind,ind)
zs3s = zeros(ind,ind)
ys3x = zeros(ind,ind)
zs3x = zeros(ind,ind)
xl₁ = zeros(ind)
yl₁ = zeros(ind)
zl = zeros(ind)
xl₂ = zeros(ind)
yl₂ = zeros(ind)
xl₃ = zeros(ind)
yl₃ = zeros(ind)
xl₄ = zeros(ind)
yl₄ = zeros(ind)
xl₅ = zeros(ind)
yl₅ = zeros(ind)
xl₆ = zeros(ind)
yl₆ = zeros(ind)
xl₇ = zeros(ind)
yl₇ = zeros(ind)
xl₈ = zeros(ind)
yl₈ = zeros(ind)
xl₉ = zeros(ind)
yl₉ = zeros(ind)
xl₁₀ = zeros(ind)
yl₁₀ = zeros(ind)
yl1ₛ = zeros(ind)
yl1ₓ = zeros(ind)
zl1ₛ = zeros(ind)
zl1ₓ = zeros(ind)
yl2ₛ = zeros(ind)
yl2ₓ = zeros(ind)
zl2ₛ = zeros(ind)
zl2ₓ = zeros(ind)
yl3ₛ = zeros(ind)
yl3ₓ = zeros(ind)
zl3ₛ = zeros(ind)
zl3ₓ = zeros(ind)
yl4ₛ = zeros(ind)
yl4ₓ = zeros(ind)
zl4ₛ = zeros(ind)
zl4ₓ = zeros(ind)
yl5ₛ = zeros(ind)
yl5ₓ = zeros(ind)
zl5ₛ = zeros(ind)
zl5ₓ = zeros(ind)
yl6ₛ = zeros(ind)
yl6ₓ = zeros(ind)
zl6ₛ = zeros(ind)
zl6ₓ = zeros(ind)
yl7ₛ = zeros(ind)
yl7ₓ = zeros(ind)
zl7ₛ = zeros(ind)
zl7ₓ = zeros(ind)
yl8ₛ = zeros(ind)
yl8ₓ = zeros(ind)
zl8ₛ = zeros(ind)
zl8ₓ = zeros(ind)
yl9ₛ = zeros(ind)
yl9ₓ = zeros(ind)
zl9ₛ = zeros(ind)
zl9ₓ = zeros(ind)
yl10ₛ = zeros(ind)
yl10ₓ = zeros(ind)
zl10ₛ = zeros(ind)
zl10ₓ = zeros(ind)
xld1 = zeros(ind)
yld1 = zeros(ind)
zld1 = zeros(ind)
xld2 = zeros(ind)
yld2 = zeros(ind)
zld2 = zeros(ind)
xld3 = zeros(ind)
yld3 = zeros(ind)
zld3 = zeros(ind)
xld4 = zeros(ind)
yld4 = zeros(ind)
zld4 = zeros(ind)
xld5 = zeros(ind)
yld5 = zeros(ind)
zld5 = zeros(ind)
xld6 = zeros(ind)
yld6 = zeros(ind)
zld6 = zeros(ind)
xld7 = zeros(ind)
yld7 = zeros(ind)
zld7 = zeros(ind)
xld8 = zeros(ind)
yld8 = zeros(ind)
zld8 = zeros(ind)
xld9 = zeros(ind)
yld9 = zeros(ind)
zld9 = zeros(ind)
xld10 = zeros(ind)
yld10 = zeros(ind)
zld10 = zeros(ind)

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
           xs1[I,J] = ξ¹
           ys1[I,J] = ξ²
           zs1[I,J] = α*u₃
           cs1[I,J] = u₃
           ys1s[I,J] = ξ²+α*u₃/(360-ξ²)*h
           zs1s[I,J] = α*u₃+5
           ys1x[I,J] = ξ²-α*u₃/(360-ξ²)*h
           zs1x[I,J] = α*u₃-5
           xl₂[J] = 0
           yl₂[J] = ξ²
           xl₃[J] = ξ¹
           yl₃[J] = ξ²
           yl2ₛ[J] = ys1s[1,J]
           zl2ₛ[J] = zs1s[1,J]
           yl2ₓ[J] = ys1x[1,J]
           zl2ₓ[J] = zs1x[1,J]
           yl3ₛ[J] = ys1s[ind,J]
           zl3ₛ[J] = zs1s[ind,J]
           yl3ₓ[J] = ys1x[ind,J]
           zl3ₓ[J] = zs1x[ind,J]
    end
    xl₁[I] = ξ¹
    yl₁[I] = 0
    zl[I] = 0
    yl1ₛ[I] = ys1s[I,1]
    zl1ₛ[I] = zs1s[I,1]
    yl1ₓ[I] = ys1x[I,1]
    zl1ₓ[I] = zs1x[I,1]
    for (i,x) in enumerate(LinRange(0.0,30, ind))
        xl₄[i]=x
        yl₄[i]=30
        xl₅[i]=x+90
        yl₅[i]=30
    end
    yl4ₛ[I] = ys1s[I,ind]
    zl4ₛ[I] = zs1s[I,ind]
    yl4ₓ[I] = ys1x[I,ind]
    zl4ₓ[I] = zs1x[I,ind]
    yl5ₛ[I] = ys1s[I,ind]
    zl5ₛ[I] = zs1s[I,ind]
    yl5ₓ[I] = ys1x[I,ind]
    zl5ₓ[I] = zs1x[I,ind]
    xld1[I] = 0
    xld2[I] = 120
    xld3[I] = 0
    xld4[I] = 30
    xld5[I] = 90
    xld6[I] = 120
    for (i,x) in enumerate(LinRange(yl1ₛ[1],yl1ₓ[1], ind))
        yld1[i] = x
    end
    for (i,x) in enumerate(LinRange(zl1ₛ[1],zl1ₓ[1,1], ind))
        zld1[i] = x
    end
    for (i,x) in enumerate(LinRange(yl1ₛ[ind],yl1ₓ[ind], ind))
        yld2[i] = x
    end
    for (i,x) in enumerate(LinRange(zl1ₛ[ind],zl1ₓ[ind], ind))
        zld2[i] = x
    end
    for (i,x) in enumerate(LinRange(yl4ₛ[1],yl4ₓ[1], ind))
        yld3[i] = x
    end
    for (i,x) in enumerate(LinRange(zl4ₛ[1],zl4ₓ[1,1], ind))
        zld3[i] = x
    end
    for (i,x) in enumerate(LinRange(yl4ₛ[ind],yl4ₓ[ind], ind))
        yld4[i] = x
    end
    for (i,x) in enumerate(LinRange(zl4ₛ[ind],zl4ₓ[ind], ind))
        zld4[i] = x
    end
    for (i,x) in enumerate(LinRange(yl5ₛ[1],yl5ₓ[1], ind))
        yld5[i] = x
    end
    for (i,x) in enumerate(LinRange(zl5ₛ[1],zl5ₓ[1,1], ind))
        zld5[i] = x
    end
    for (i,x) in enumerate(LinRange(yl5ₛ[ind],yl5ₓ[ind], ind))
        yld6[i] = x
    end
    for (i,x) in enumerate(LinRange(zl5ₛ[ind],zl5ₓ[ind], ind))
        zld6[i] = x
    end
end

for i in 1:ind
    for j in 1:ind
        Δy=50/(ind-1)
        Δx=(0.3*(i-1)Δy+15)/(ind-1)*2
        x₂ = 60 + (j-(ind+1)/2)*Δx
        y₂ = 80 - (i-1)*Δy
        indices = sp(x₂,y₂,0.0)
        N = zeros(length(indices))
        data = Dict([:x=>(2,[x₂]),:y=>(2,[y₂]),:z=>(2,[0.0]),:𝝭=>(4,N),:𝗠=>(0,𝗠)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set𝝭!(ap)
        u₃ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
        end
        xs2[i,j] = x₂
        ys2[i,j] = y₂
        zs2[i,j] = α*u₃
        cs2[i,j] = u₃
        ys2s[i,j] = y₂+α*u₃/(360-y₂)*h
        zs2s[i,j] = α*u₃+5
        ys2x[i,j] = y₂-α*u₃/(360-y₂)*h
        zs2x[i,j] = α*u₃-5
        xl₆[i]=60-(0.3*(i-1)Δy+15)
        yl₆[i]=y₂
        xl₇[i]=60+(0.3*(i-1)Δy+15)
        yl₇[i]=y₂
    end
    yl6ₛ[i] = ys2s[i,1]
    zl6ₛ[i] = zs2s[i,1]
    yl6ₓ[i] = ys2x[i,1]
    zl6ₓ[i] = zs2x[i,1]   
    yl7ₛ[i] = ys2s[i,ind]
    zl7ₛ[i] = zs2s[i,ind]
    yl7ₓ[i] = ys2x[i,ind]
    zl7ₓ[i] = zs2x[i,ind]
for i in 1:ind
    for j in 1:ind
        Δy=280/(ind-1)
        Δx=(-0.2*(i-1)Δy+71)/(ind-1)*2
        x₃ = 60 + (j-(ind+1)/2)*Δx
        y₃ = 360 - (i-1)*Δy
        indices = sp(x₃,y₃,0.0)
        N = zeros(length(indices))
        data = Dict([:x=>(2,[x₃]),:y=>(2,[y₃]),:z=>(2,[0.0]),:𝝭=>(4,N),:𝗠=>(0,𝗠)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set𝝭!(ap)
        u₃ = 0.0
        for (i,xᵢ) in enumerate(𝓒)
            u₃ += N[i]*xᵢ.d
        end
        xs3[i,j] = x₃
        ys3[i,j] = y₃
        zs3[i,j] = α*u₃
        cs3[i,j] = u₃
        ys3s[i,j] = y₃+α*u₃/(360.1-y₃)*h
        zs3s[i,j] = α*u₃+5
        ys3x[i,j] = y₃-α*u₃/(360.1-y₃)*h
        zs3x[i,j] = α*u₃-5
        xl₈[i]=60+(-0.2*(i-1)Δy+71)
        yl₈[i]=y₃
        xl₉[i]=60-(-0.2*(i-1)Δy+71)
        yl₉[i]=y₃
        yl10ₛ[j] = ys3s[1,j]
        zl10ₛ[j] = zs3s[1,j]
        yl10ₓ[j] = ys3x[1,j]
        zl10ₓ[j] = zs3x[1,j]   
    end
    for (i,x) in enumerate(LinRange(-11,131, ind))
        xl₁₀[i]=x
        yl₁₀[i]=360
    end
    yl8ₛ[i] = ys3s[i,1]
    zl8ₛ[i] = zs3s[i,1]
    yl8ₓ[i] = ys3x[i,1]
    zl8ₓ[i] = zs3x[i,1]   
    yl9ₛ[i] = ys3s[i,ind]
    zl9ₛ[i] = zs3s[i,ind]
    yl9ₓ[i] = ys3x[i,ind]
    zl9ₓ[i] = zs3x[i,ind]
    xld7[i] = -11
    xld8[i] = 45
    xld9[i] = 131
    xld10[i] = 75
    for (I,x) in enumerate(LinRange(yl8ₛ[1],yl8ₓ[1], ind))
        yld7[I] = x
    end
    for (I,x) in enumerate(LinRange(zl8ₛ[1],zl8ₓ[1,1], ind))
        zld7[I] = x
    end
    for (I,x) in enumerate(LinRange(yl8ₛ[ind],yl8ₓ[ind], ind))
        yld8[I] = x
    end
    for (I,x) in enumerate(LinRange(zl8ₛ[ind],zl8ₓ[ind], ind))
        zld8[I] = x
    end
    for (I,x) in enumerate(LinRange(yl9ₛ[1],yl9ₓ[1], ind))
        yld9[I] = x
    end
    for (I,x) in enumerate(LinRange(zl9ₛ[1],zl9ₓ[1,1], ind))
        zld9[I] = x
    end
    for (I,x) in enumerate(LinRange(yl9ₛ[ind],yl9ₓ[ind], ind))
        yld10[I] = x
    end
    for (I,x) in enumerate(LinRange(zl9ₛ[ind],zl9ₓ[ind], ind))
        zld10[I] = x
    end
end
end

fig = Figure()
ax = Axis3(fig[1, 1], aspect = :data, azimuth = 0.75*pi, elevation = 0.80*pi)

hidespines!(ax)
hidedecorations!(ax)
s = surface!(ax,xs1,ys1,zs1, color=cs1, colormap=:redsblues, colorrange = (-0.11,0))
# s = surface!(ax,xs1s,ys1s,zs1s, color=cs1, colormap=:redsblues, colorrange = (-0.11,0))
# s = surface!(ax,xs1x,ys1x,zs1x, color=cs1, colormap=:redsblues, colorrange = (-0.11,0))
s = surface!(ax,xs2,ys2,zs2, color=cs2, colormap=:redsblues, colorrange = (-0.11,0))
# s = surface!(ax,xs2s,ys2s,zs2s, color=cs2, colormap=:redsblues, colorrange = (-0.11,0))
# s = surface!(ax,xs2x,ys2x,zs2x, color=cs2, colormap=:redsblues, colorrange = (-0.11,0))
s = surface!(ax,xs3,ys3,zs3, color=cs3, colormap=:redsblues, colorrange = (-0.11,0))
# s = surface!(ax,xs3s,ys3s,zs3s, color=cs3, colormap=:redsblues, colorrange = (-0.11,0))
# s = surface!(ax,xs3x,ys3x,zs3x, color=cs3, colormap=:redsblues, colorrange = (-0.11,0))
lines!(ax,xl₁,yl₁,zl,color=:purple)
lines!(ax,xl₂,yl₂,zl,color=:purple)
lines!(ax,xl₃,yl₃,zl,color=:purple)
lines!(ax,xl₄,yl₄,zl,color=:purple)
lines!(ax,xl₅,yl₅,zl,color=:purple)
lines!(ax,xl₆,yl₆,zl,color=:purple)
lines!(ax,xl₇,yl₇,zl,color=:purple)
lines!(ax,xl₈,yl₈,zl,color=:purple)
lines!(ax,xl₉,yl₉,zl,color=:purple)
lines!(ax,xl₁₀,yl₁₀,zl,color=:purple)
lines!(ax,xl₁,yl1ₛ,zl1ₛ,color=:gray)
lines!(ax,xl₁,yl1ₓ,zl1ₓ,color=:gray)
lines!(ax,xl₂,yl2ₛ,zl2ₛ,color=:gray)
lines!(ax,xl₂,yl2ₓ,zl2ₓ,color=:gray)
lines!(ax,xl₃,yl3ₛ,zl3ₛ,color=:gray)
lines!(ax,xl₃,yl3ₓ,zl3ₓ,color=:gray)
lines!(ax,xl₄,yl4ₛ,zl4ₛ,color=:gray)
lines!(ax,xl₄,yl4ₓ,zl4ₓ,color=:gray)
lines!(ax,xl₅,yl5ₛ,zl5ₛ,color=:gray)
lines!(ax,xl₅,yl5ₓ,zl5ₓ,color=:gray)
lines!(ax,xl₆,yl6ₛ,zl6ₛ,color=:gray)
lines!(ax,xl₆,yl6ₓ,zl6ₓ,color=:gray)
lines!(ax,xl₇,yl7ₛ,zl7ₛ,color=:gray)
lines!(ax,xl₇,yl7ₓ,zl7ₓ,color=:gray)
lines!(ax,xl₈,yl8ₛ,zl8ₛ,color=:gray)
lines!(ax,xl₈,yl8ₓ,zl8ₓ,color=:gray)
lines!(ax,xl₉,yl9ₛ,zl9ₛ,color=:gray)
lines!(ax,xl₉,yl9ₓ,zl9ₓ,color=:gray)
lines!(ax,xl₁₀,yl10ₛ,zl10ₛ,color=:gray)
lines!(ax,xl₁₀,yl10ₓ,zl10ₓ,color=:gray)
lines!(ax,xld1,yld1,zld1,color=:gray)
lines!(ax,xld2,yld2,zld2,color=:gray)
lines!(ax,xld3,yld3,zld3,color=:gray)
lines!(ax,xld4,yld4,zld4,color=:gray)
lines!(ax,xld5,yld5,zld5,color=:gray)
lines!(ax,xld6,yld6,zld6,color=:gray)
lines!(ax,xld7,yld7,zld7,color=:gray)
lines!(ax,xld8,yld8,zld8,color=:gray)
lines!(ax,xld9,yld9,zld9,color=:gray)
lines!(ax,xld10,yld10,zld10,color=:gray)

Colorbar(fig[2, 1], s, vertical = false)

fig
