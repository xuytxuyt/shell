using ApproxOperator, JLD, GLMakie, Tensors

import Gmsh: gmsh
import BenchmarkExample: BenchmarkExample

𝑅 = BenchmarkExample.PatchTestThinShell.𝑅
𝐿 = BenchmarkExample.PatchTestThinShell.𝐿
𝜃 = BenchmarkExample.PatchTestThinShell.𝜃
E = BenchmarkExample.PatchTestThinShell.𝐸
ν = BenchmarkExample.PatchTestThinShell.𝜈
h = BenchmarkExample.PatchTestThinShell.ℎ
Dᵇ = E*h^3/12/(1-ν^2)

n = 2
u(x) = Tensors.Vec{3}(((x[1]+2*x[2])^n,(3*x[1]+4*x[2])^n,(5*x[1]+6*x[2])^n))
css = BenchmarkExample.cylindricalCoordinate(𝑅)
vss = BenchmarkExample.PatchTestThinShell.variables(css,u)

## import nodes
gmsh.initialize()
gmsh.open("msh/patchtest.msh")
entities = getPhysicalGroups()
nodes = get𝑿ᵢ()
x = nodes.x
y = nodes.y
z = nodes.z
sp = RegularGrid(x,y,z,n = 1,γ = 5)
nₚ = length(nodes)
s = 2.5*0.1*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
elements = Dict{String,Vector{ApproxOperator.AbstractElement}}()
elements["Ω"] = getElements(nodes,entities["Ω"])
elements["Γ¹"] = getElements(nodes,entities["Γ¹"])
elements["Γ²"] = getElements(nodes,entities["Γ²"])
elements["Γ³"] = getElements(nodes,entities["Γ³"])
elements["Γ⁴"] = getElements(nodes,entities["Γ⁴"])
elements["∂Ω"] = elements["Γ¹"]∪elements["Γ²"]∪elements["Γ³"]∪elements["Γ⁴"]
gmsh.finalize()

type = ReproducingKernel{:Quadratic2D,:□,:QuinticSpline}
# 𝗠 = zeros(21)
# type = ReproducingKernel{:Cubic2D,:□,:CubicSpline}
𝗠 = zeros(55)

# filename = "patchtest_gauss_penalty"
# filename = "patchtest_gauss_nitsche"
# filename = "patchtest_mix_penalty"
# filename = "patchtest_mix_nitsche"
filename = "patchtest_mix_hr"
ds = Dict(load("jld/"*filename*".jld"))
push!(nodes,:d₁=>ds["d₁"],:d₂=>ds["d₂"],:d₃=>ds["d₃"])

ind = 21
xs = zeros(ind,ind)
ys = zeros(ind,ind)
zs = zeros(ind,ind)
cs = zeros(ind,ind)
us = zeros(ind,ind)
vs = zeros(ind,ind)
ws = zeros(ind,ind)
for (I,ξ¹) in enumerate(LinRange(0.0, 𝜃*𝑅, ind))
    for (J,ξ²) in enumerate(LinRange(0.0, 𝐿, ind))
        indices = sp(ξ¹,ξ²,0.0)
        N = zeros(length(indices))
        B₁ = zeros(length(indices))
        B₂ = zeros(length(indices))
        B₁₁ = zeros(length(indices))
        B₁₂ = zeros(length(indices))
        B₂₂ = zeros(length(indices))
        B₁₁₁ = zeros(length(indices))
        B₁₁₂ = zeros(length(indices))
        B₁₂₂ = zeros(length(indices))
        B₂₂₂ = zeros(length(indices))
        data = Dict([:x=>(2,[ξ¹]),:y=>(2,[ξ²]),:z=>(2,[0.0]),:𝝭=>(4,N),:∂𝝭∂x=>(4,B₁),:∂𝝭∂y=>(4,B₂),:∂²𝝭∂x²=>(4,B₁₁),:∂²𝝭∂x∂y=>(4,B₁₂),:∂²𝝭∂y²=>(4,B₂₂),:∂³𝝭∂x³=>(4,B₁₁₁),:∂³𝝭∂x²∂y=>(4,B₁₁₂),:∂³𝝭∂x∂y²=>(4,B₁₂₂),:∂³𝝭∂y³=>(4,B₂₂₂),:𝗠=>(0,𝗠)])
        𝓒 = [nodes[k] for k in indices]
        𝓖 = [𝑿ₛ((𝑔=1,𝐺=1,𝐶=1,𝑠=0),data)]
        ap = type(𝓒,𝓖)
        set∇̂³𝝭!(ap)
        x = Tensors.Vec{3}((ξ¹,ξ²,0.0))
        u₁ = 0.0
        u₂ = 0.0
        u₃ = 0.0
        𝑴ex = vss.𝑴(x)
        𝑴 = zeros(3)
        a¹¹ = css.a¹¹(x)
        a¹² = css.a¹²(x)
        a²² = css.a²²(x)
        𝒂₁₍₁₎ = css.𝒂₁(x)[1]
        𝒂₁₍₂₎ = css.𝒂₁(x)[2]
        𝒂₁₍₃₎ = css.𝒂₁(x)[3]
        𝒂₂₍₁₎ = css.𝒂₂(x)[1]
        𝒂₂₍₂₎ = css.𝒂₂(x)[2]
        𝒂₂₍₃₎ = css.𝒂₂(x)[3]
        𝒂₃₍₁₎ = css.𝒂₃(x)[1]
        𝒂₃₍₂₎ = css.𝒂₃(x)[2]
        𝒂₃₍₃₎ = css.𝒂₃(x)[3]
        Γ¹₁₁ = css.Γ¹₁₁(x)
        Γ¹₂₂ = css.Γ¹₂₂(x)
        Γ¹₁₂ = css.Γ¹₁₂(x)
        Γ²₁₁ = css.Γ²₁₁(x)
        Γ²₂₂ = css.Γ²₂₂(x)
        Γ²₁₂ = css.Γ²₁₂(x)
        𝐃 = Tensors.@SArray [                  a¹¹*a¹¹ ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                            a¹¹*a¹²;
                     ν*a¹¹*a²² + (1-ν)*a¹²*a¹²                   a²²*a²²                            a²²*a¹²;
                                       a¹¹*a¹²                   a²²*a¹² 0.5*((1-ν)*a¹¹*a²² + (1+ν)*a¹²*a¹²)]
        for (i,xᵢ) in enumerate(𝓒)
            u₁ += N[i]*xᵢ.d₁
            u₂ += N[i]*xᵢ.d₂
            u₃ += N[i]*xᵢ.d₃
            Bᵇ₁₁ = Γ¹₁₁*B₁[i]+Γ²₁₁*B₂[i]-B₁₁[i]
            Bᵇ₂₂ = Γ¹₂₂*B₁[i]+Γ²₂₂*B₂[i]-B₂₂[i]
            Bᵇ₁₂ = Γ¹₁₂*B₁[i]+Γ²₁₂*B₂[i]-B₁₂[i]
            𝐁ᵇᵢ = Tensors.@SArray [  Bᵇ₁₁*𝒂₃₍₁₎   Bᵇ₁₁*𝒂₃₍₂₎   Bᵇ₁₁*𝒂₃₍₃₎;
                             Bᵇ₂₂*𝒂₃₍₁₎   Bᵇ₂₂*𝒂₃₍₂₎   Bᵇ₂₂*𝒂₃₍₃₎;
                           2*Bᵇ₁₂*𝒂₃₍₁₎ 2*Bᵇ₁₂*𝒂₃₍₂₎ 2*Bᵇ₁₂*𝒂₃₍₃₎]
            dᵢ = Tensors.Vec{3}((xᵢ.d₁,xᵢ.d₂,xᵢ.d₃))
            𝑴 .+= Dᵇ*𝐃*𝐁ᵇᵢ*dᵢ
        end
        us[I,J] = u₁
        vs[I,J] = u₂
        ws[I,J] = u₃
        xs[I,J] = 𝑅*sin(ξ¹/𝑅)
        ys[I,J] = ξ²
        zs[I,J] = 𝑅*cos(ξ¹/𝑅)
        cs[I,J] = 𝑴ex[3]
        # cs[I,J] = 𝑴[3]
    end
end
# colorrange1 = (0.0,9.0)
# colorrange2 = (0.0,49.0)
# colorrange3 = (0.0,121.0)
# colorrange1 = (-0.005,-0.002)
# colorrange2 = (-0.0055,-0.0025)
colorrange3 = (-0.005,-0.003)

fig = Figure()
ax = Axis3(fig[1, 1], perspectiveness = 0.8, aspect = :data, azimuth = -0.25*pi, elevation = 0.2*pi, xlabel = " ", ylabel = " ", zlabel = " ", xticksvisible = false,xticklabelsvisible=false, yticksvisible = false, yticklabelsvisible=false, zticksvisible = false, zticklabelsvisible=false, protrusions = (0.,0.,0.,0.))
Colorbar(fig[1,2], colormap=:lightrainbow, size = 16, height = 250, tellheight=false, ticklabelsvisible=false, ticksvisible=false, minorticksvisible = false, limits=colorrange3)

hidespines!(ax)
# hidedecorations!(ax)
# s = surface!(ax,xs,ys,zs, color=cs, colormap=:lightrainbow, colorrange = colorrange3)
s = surface!(ax,xs,ys,zs, color=cs, colormap=:lightrainbow)

for elm in elements["∂Ω"]
    ξ¹ = [x.x for x in elm.𝓒]
    ξ² = [x.y for x in elm.𝓒]
    x = [𝑅*sin(ξ/𝑅) for ξ in ξ¹]
    y = ξ²
    z = [𝑅*cos(ξ/𝑅) for ξ in ξ¹]
    lines!(x,y,z,linewidth = 1.5, color = :black)
end

# save("./png/"*filename*".png",fig)
save("./png/patchtest_exact.png",fig)
fig
