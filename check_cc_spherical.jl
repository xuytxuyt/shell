using Revise
using ApproxOperator, Tensors, BenchmarkExample, LinearAlgebra

include("import_prescrible_ops.jl")
include("import_Spherical_shell.jl")
ndiv = 8
elements, nodes = import_spherical_mix("msh/sphericalshell_"*string(ndiv)*".msh");

nₘ = 21
𝑅 = BenchmarkExample.SphericalShell.𝑅
E = BenchmarkExample.SphericalShell.𝐸
ν = BenchmarkExample.SphericalShell.𝜈
h = BenchmarkExample.SphericalShell.ℎ
𝜃 =  BenchmarkExample.SphericalShell.𝜃₂
𝐹 = BenchmarkExample.SphericalShell.𝐹

cs = BenchmarkExample.sphericalCoordinate(𝑅)
nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])*3)
s = 2.5*𝑅*𝜃/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

eval(prescribeForMix)
eval(prescribeVariables)

set𝝭!(elements["Ω"])
set∇²𝝭!(elements["Ωₚ"])
set∇𝝭!(elements["Γₚ"])
set∇𝝭!(elements["Ωₘ"])
set∇𝝭!(elements["Γₘ"])
set∇𝝭!(elements["Γʳ"])
set∇𝝭!(elements["Γˡ"])
set𝝭!(elements["𝐴"])
set𝝭!(elements["𝐵"])

n = 0
uex(x) = Vec{3}(((x[1]+2*x[2])^n,(3*x[1]+4*x[2])^n,(5*x[1]+6*x[2])^n))
∂uex(x) = gradient(uex,x)
∂²uex(x) = gradient(∂uex,x)
dᵛ = zeros(3*nₚ)
d₁ = zeros(nₚ)
d₂ = zeros(nₚ)
d₃ = zeros(nₚ)
for (I,node) in enumerate(nodes)
    x = Vec{3}((node.x,node.y,node.z))
    u = uex(x)
    dᵛ[3*I-2] = u[1]
    dᵛ[3*I-1] = u[2]
    dᵛ[3*I]   = u[3]
    d₁[I] = u[1]
    d₂[I] = u[2]
    d₃[I] = u[3]
end
push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)

eval(opsMix)

kᵋᵋ = zeros(6*nᵥ,6*nᵥ)
kᴺᵋ = zeros(6*nᵥ,6*nᵥ)
kᴺᵛ = zeros(6*nᵥ,3*nₚ)
kᵏᵏ = zeros(9*nᵥ,9*nᵥ)
kᴹᵏ = zeros(9*nᵥ,9*nᵥ)
kᴹᵛ = zeros(9*nᵥ,3*nₚ)

ops[1](elements["Ω"],kᵋᵋ)
ops[2](elements["Ω"],kᴺᵋ)
ops[3](elements["Γₚ"],elements["Γₘ"],kᴺᵛ)
ops[4](elements["Ωₚ"],elements["Ωₘ"],kᴺᵛ)
# ops[3](elements["Γₚ"][1:3],elements["Γₘ"][1:3],kᴺᵛ)
# ops[4](elements["Ωₚ"][1:1],elements["Ωₘ"][1:1],kᴺᵛ)

ops[5](elements["Ω"],kᵏᵏ)
ops[6](elements["Ω"],kᴹᵏ)
ops[7](elements["Γₚ"],elements["Γₘ"],kᴹᵛ)
ops[8](elements["Γₚ"],elements["Γₘ"],kᴹᵛ)
ops[9](elements["Γₚ"],elements["Γₘ"],kᴹᵛ)
ops[10](elements["Ωₚ"],elements["Ωₘ"],kᴹᵛ)
# ops[7](elements["Γₚ"][1:3],elements["Γₘ"][1:3],kᴹᵛ)
# ops[8](elements["Γₚ"][1:3],elements["Γₘ"][1:3],kᴹᵛ)
# ops[9](elements["Γₚ"][1:3],elements["Γₘ"][1:3],kᴹᵛ)
# ops[10](elements["Ωₚ"][1:1],elements["Ωₘ"][1:1],kᴹᵛ)

kᵇ = - kᴺᵋ\kᴺᵛ
kᵐ = - kᴹᵏ\kᴹᵛ

fᵇ = kᵇ*dᵛ
fᵐ = kᵐ*dᵛ

err = 0.0
for a in elements["Ωₚ"]
    𝓒 = a.𝓒
    𝓖 = a.𝓖
    for ξ in 𝓖
        Γ¹₁₁ = ξ.Γ¹₁₁
        Γ²₁₁ = ξ.Γ²₁₁
        Γ¹₂₂ = ξ.Γ¹₂₂
        Γ²₂₂ = ξ.Γ²₂₂
        Γ¹₁₂ = ξ.Γ¹₁₂
        Γ²₁₂ = ξ.Γ²₁₂
        N = ξ[:𝝭]
        x = Vec{3}((ξ.x,ξ.y,ξ.z))
        ∂₁u = ∂uex(x)[:,1]
        ∂₂u = ∂uex(x)[:,2]
        ∂₁₁u = Γ¹₁₁*∂uex(x)[:,1] + Γ²₁₁*∂uex(x)[:,2] - ∂²uex(x)[:,1,1]
        ∂₂₂u = Γ¹₂₂*∂uex(x)[:,1] + Γ²₂₂*∂uex(x)[:,2] - ∂²uex(x)[:,2,2]
        ∂₁₂u = 2*(Γ¹₁₂*∂uex(x)[:,1] + Γ²₁₂*∂uex(x)[:,2] - ∂²uex(x)[:,1,2])
        ∂₁uʰ = zeros(3)
        ∂₂uʰ = zeros(3)
        ∂₁₁uʰ = zeros(3)
        ∂₂₂uʰ = zeros(3)
        ∂₁₂uʰ = zeros(3)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            ∂₁uʰ[1] += N[i]*fᵇ[6*I-5]
            ∂₁uʰ[2] += N[i]*fᵇ[6*I-4]
            ∂₁uʰ[3] += N[i]*fᵇ[6*I-4]
            ∂₂uʰ[1] += N[i]*fᵇ[6*I-2]
            ∂₂uʰ[2] += N[i]*fᵇ[6*I-1]
            ∂₂uʰ[3] += N[i]*fᵇ[6*I]

            ∂₁₁uʰ[1] += N[i]*fᵐ[9*I-8]
            ∂₁₁uʰ[2] += N[i]*fᵐ[9*I-7]
            ∂₁₁uʰ[3] += N[i]*fᵐ[9*I-6]
            ∂₂₂uʰ[1] += N[i]*fᵐ[9*I-5]
            ∂₂₂uʰ[2] += N[i]*fᵐ[9*I-4]
            ∂₂₂uʰ[3] += N[i]*fᵐ[9*I-3]
            ∂₁₂uʰ[1] += N[i]*fᵐ[9*I-2]
            ∂₁₂uʰ[2] += N[i]*fᵐ[9*I-1]
            ∂₁₂uʰ[3] += N[i]*fᵐ[9*I]
        end
        # println(∂₁u[3])
        # println(∂₂u[3])
        # println(∂₂uʰ[1])
        # println(∂₁uʰ[1])
        # global err += abs(∂₁u[3]-∂₁uʰ[3])
        # global err += abs(∂₂u[2]-∂₂uʰ[2])
        # global err += abs(∂₁₁u[3]-∂₁₁uʰ[3])
        # global err += abs(∂₂₂u[3]-∂₂₂uʰ[3])
        global err += abs(∂₁₂u[2]-∂₁₂uʰ[2])
    end
end
println(err)

# f1 = zeros(3)
# f2 = zeros(3)
# kt = zeros(3,3)
# f = 0.0
# for a in elements["Γₚ"][1:3]
#     𝓒 = a.𝓒
#     𝓖 = a.𝓖
#     for ξ in 𝓖
#         𝑤 = ξ.𝑤
#         n₁ = ξ.n₁
#         n₂ = ξ.n₂
#         p = [1.0,ξ.x,ξ.y]
#         f1 .+= - p*n₁*𝑤
#         f2 .+= - p*n₂*𝑤
#         # global f -= ξ.x*n₁*𝑤
#     end
# end
# for a in elements["Ω"][1:1]
#     𝓒 = a.𝓒
#     𝓖 = a.𝓖
#     for ξ in 𝓖
#         𝑤 = ξ.𝑤
#         Γ¹₁₁ = ξ.Γ¹₁₁
#         Γ¹₁₂ = ξ.Γ¹₁₂
#         Γ¹₂₂ = ξ.Γ¹₂₂
#         Γ²₁₁ = ξ.Γ²₁₁
#         Γ²₁₂ = ξ.Γ²₁₂
#         Γ²₂₂ = ξ.Γ²₂₂
#         Γᵞᵧ₁ = Γ¹₁₁ + Γ²₁₂
#         Γᵞᵧ₂ = Γ¹₁₂ + Γ²₂₂
#         p = [1.0,ξ.x,ξ.y]
#         p₁ = [0.0,1.0,0.0]
#         p₂ = [0.0,0.0,1.0]
#         kt .+= p*p'*𝑤
#         f1 .+= (p₁ + Γᵞᵧ₁*p)*𝑤
#         f2 .+= (p₂ + Γᵞᵧ₂*p)*𝑤
#         # global f += 𝑤
#         # global f -= ξ.x*(Γᵞᵧ₁)*𝑤
#         # global f += ξ.w
#         # println(𝑤)
#         println(f1)
#     end
#     # println(a.𝐴)
#     # global f += - a.𝐴
# end

# println(kt\f1)
# println(f)

# f1 = kᴺᵛ*dᵛ
# f2 = zeros(6)
# f3 = zeros(6)
f1 = kᴹᵛ*dᵛ
f2 = zeros(27)
f3 = zeros(27)
temp = 0.0
for a in elements["Γₚ"][1:3]
# for a in elements["Γₚ"]
    𝓒 = a.𝓒
    𝓖 = a.𝓖
    for ξ in 𝓖
        x_ = Vec{3}((ξ.x,ξ.y,0.0))
        𝑤 = ξ.𝑤
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        n¹ = ξ.n¹
        n² = ξ.n²
        s₁ = ξ.s₁
        s₂ = ξ.s₂
        s¹ = ξ.s¹
        s² = ξ.s²
        Δ = ξ.Δ

        ∂₁n₁ = ξ.∂₁n₁
        ∂₂n₁ = ξ.∂₂n₁
        ∂₁n₂ = ξ.∂₁n₂
        ∂₂n₂ = ξ.∂₂n₂
        ∂₁s₁ = ξ.∂₁s₁
        ∂₂s₁ = ξ.∂₂s₁
        ∂₁s₂ = ξ.∂₁s₂
        ∂₂s₂ = ξ.∂₂s₂
        Γ¹₁₁ = ξ.Γ¹₁₁
        Γ¹₁₂ = ξ.Γ¹₁₂
        Γ¹₂₂ = ξ.Γ¹₂₂
        Γ²₁₁ = ξ.Γ²₁₁
        Γ²₁₂ = ξ.Γ²₁₂
        Γ²₂₂ = ξ.Γ²₂₂
        Γᵞᵧ₁n₁ = Γ¹₁₁*n₁ + Γ²₁₂*n₁
        Γᵞᵧ₂n₂ = Γ¹₁₂*n₂ + Γ²₂₂*n₂
        Γᵞᵧ₂n₁ = Γ¹₁₂*n₁ + Γ²₂₂*n₁
        Γᵞᵧ₁n₂ = Γ¹₁₁*n₂ + Γ²₁₂*n₂
        Γᵞ₁₁nᵧ = Γ¹₁₁*n₁ + Γ²₁₁*n₂
        Γᵞ₂₂nᵧ = Γ¹₂₂*n₁ + Γ²₂₂*n₂
        Γᵞ₁₂nᵧ = Γ¹₁₂*n₁ + Γ²₁₂*n₂
        s₁n₁s¹ = s₁*n₁*s¹
        s₁n₁s² = s₁*n₁*s²
        s₂n₂s¹ = s₂*n₂*s¹
        s₂n₂s² = s₂*n₂*s²
        s₁n₂s¹ = s₁*n₂*s¹
        s₁n₂s² = s₁*n₂*s²
        s₂n₁s¹ = s₂*n₁*s¹
        s₂n₁s² = s₂*n₁*s²
        ∂ᵧs₁n₁sᵞ = ∂₁s₁*n₁*s¹ + ∂₂s₁*n₁*s²
        ∂ᵧs₂n₂sᵞ = ∂₁s₂*n₂*s¹ + ∂₂s₂*n₂*s²
        ∂ᵧs₁n₂sᵞ = ∂₁s₁*n₂*s¹ + ∂₂s₁*n₂*s²
        ∂ᵧs₂n₁sᵞ = ∂₁s₂*n₁*s¹ + ∂₂s₂*n₁*s²
        s₁∂ᵧn₁sᵞ = s₁*∂₁n₁*s¹ + s₁*∂₂n₁*s²
        s₂∂ᵧn₂sᵞ = s₂*∂₁n₂*s¹ + s₂*∂₂n₂*s²
        s₁∂ᵧn₂sᵞ = s₁*∂₁n₂*s¹ + s₁*∂₂n₂*s²
        s₂∂ᵧn₁sᵞ = s₂*∂₁n₁*s¹ + s₂*∂₂n₁*s²
        # println(n₁)
        # println(n₁*n¹+n₂*n²)
        p = [1,ξ.x,ξ.y]
        p₁ = [0.0,1.0,0.0]
        p₂ = [0.0,0.0,1.0]
        # f2[1:2:6] .-= p.*n₁.*𝑤
        # f2[2:2:6] .-= p.*n₂.*𝑤
        # f2[1:3:9] .-= ((p₁.*s¹+p₂.*s²)*s₁*n₁ + (∂₁s₁*s¹+∂₂s₁*s²).*p.*n₁ + s₁.*p*(∂₁n₁*s¹+∂₂n₁*s²) + p₁.*n₁ + Γᵞ₁₁nᵧ.*p + Γᵞᵧ₁n₁.*p)*𝑤
        f2[1:9:27] .-= ((p₁.*s¹+p₂.*s²)*s₁*n₁*𝑤 + (∂₁s₁*s¹+∂₂s₁*s²).*p.*n₁*𝑤 + s₁.*p*(∂₁n₁*s¹+∂₂n₁*s²)*𝑤 + p₁.*n₁*𝑤 + Γᵞ₁₁nᵧ.*p*𝑤 + Γᵞᵧ₁n₁.*p*𝑤)*uex(x_)[1]
        f2[2:9:27] .-= ((p₁.*s¹+p₂.*s²)*s₁*n₁*𝑤 + (∂₁s₁*s¹+∂₂s₁*s²).*p.*n₁*𝑤 + s₁.*p*(∂₁n₁*s¹+∂₂n₁*s²)*𝑤 + p₁.*n₁*𝑤 + Γᵞ₁₁nᵧ.*p*𝑤 + Γᵞᵧ₁n₁.*p*𝑤)*uex(x_)[2]
        f2[3:9:27] .-= ((p₁.*s¹+p₂.*s²)*s₁*n₁*𝑤 + (∂₁s₁*s¹+∂₂s₁*s²).*p.*n₁*𝑤 + s₁.*p*(∂₁n₁*s¹+∂₂n₁*s²)*𝑤 + p₁.*n₁*𝑤 + Γᵞ₁₁nᵧ.*p*𝑤 + Γᵞᵧ₁n₁.*p*𝑤)*uex(x_)[3]
        f2[4:9:27] .-= ((p₁.*s¹+p₂.*s²)*s₂*n₂*𝑤 + (∂₁s₂*s¹+∂₂s₂*s²).*p.*n₂*𝑤 + s₂.*p*(∂₁n₂*s¹+∂₂n₂*s²)*𝑤 + p₂.*n₂*𝑤 + Γᵞ₂₂nᵧ.*p*𝑤 + Γᵞᵧ₂n₂.*p*𝑤)*uex(x_)[1]
        f2[5:9:27] .-= ((p₁.*s¹+p₂.*s²)*s₂*n₂*𝑤 + (∂₁s₂*s¹+∂₂s₂*s²).*p.*n₂*𝑤 + s₂.*p*(∂₁n₂*s¹+∂₂n₂*s²)*𝑤 + p₂.*n₂*𝑤 + Γᵞ₂₂nᵧ.*p*𝑤 + Γᵞᵧ₂n₂.*p*𝑤)*uex(x_)[2]
        f2[6:9:27] .-= ((p₁.*s¹+p₂.*s²)*s₂*n₂*𝑤 + (∂₁s₂*s¹+∂₂s₂*s²).*p.*n₂*𝑤 + s₂.*p*(∂₁n₂*s¹+∂₂n₂*s²)*𝑤 + p₂.*n₂*𝑤 + Γᵞ₂₂nᵧ.*p*𝑤 + Γᵞᵧ₂n₂.*p*𝑤)*uex(x_)[3]
        f2[7:9:27] .-= ((p₁.*s¹+p₂.*s²)*(s₁*n₂ + s₂*n₁)*𝑤 + (∂₁s₁*s¹+∂₂s₁*s²).*p.*n₂*𝑤 + (∂₁s₂*s¹+∂₂s₂*s²).*p.*n₁*𝑤 + s₁.*p*(∂₁n₂*s¹+∂₂n₂*s²)*𝑤 + s₂.*p*(∂₁n₁*s¹+∂₂n₁*s²)*𝑤 + (p₁.*n₂ + p₂.*n₁)*𝑤 + 2*Γᵞ₁₂nᵧ.*p*𝑤 + (Γᵞᵧ₁n₂+Γᵞᵧ₂n₁).*p*𝑤)*uex(x_)[1]
        f2[8:9:27] .-= ((p₁.*s¹+p₂.*s²)*(s₁*n₂ + s₂*n₁)*𝑤 + (∂₁s₁*s¹+∂₂s₁*s²).*p.*n₂*𝑤 + (∂₁s₂*s¹+∂₂s₂*s²).*p.*n₁*𝑤 + s₁.*p*(∂₁n₂*s¹+∂₂n₂*s²)*𝑤 + s₂.*p*(∂₁n₁*s¹+∂₂n₁*s²)*𝑤 + (p₁.*n₂ + p₂.*n₁)*𝑤 + 2*Γᵞ₁₂nᵧ.*p*𝑤 + (Γᵞᵧ₁n₂+Γᵞᵧ₂n₁).*p*𝑤)*uex(x_)[2]
        f2[9:9:27] .-= ((p₁.*s¹+p₂.*s²)*(s₁*n₂ + s₂*n₁)*𝑤 + (∂₁s₁*s¹+∂₂s₁*s²).*p.*n₂*𝑤 + (∂₁s₂*s¹+∂₂s₂*s²).*p.*n₁*𝑤 + s₁.*p*(∂₁n₂*s¹+∂₂n₂*s²)*𝑤 + s₂.*p*(∂₁n₁*s¹+∂₂n₁*s²)*𝑤 + (p₁.*n₂ + p₂.*n₁)*𝑤 + 2*Γᵞ₁₂nᵧ.*p*𝑤 + (Γᵞᵧ₁n₂+Γᵞᵧ₂n₁).*p*𝑤)*uex(x_)[3]
        f2[1:9:27] .+= p*n₁*n₁ * (∂uex(x_)[1,1]*n¹ + ∂uex(x_)[1,2]*n²)*𝑤
        f2[2:9:27] .+= p*n₁*n₁ * (∂uex(x_)[2,1]*n¹ + ∂uex(x_)[2,2]*n²)*𝑤
        f2[3:9:27] .+= p*n₁*n₁ * (∂uex(x_)[3,1]*n¹ + ∂uex(x_)[3,2]*n²)*𝑤
        f2[4:9:27] .+= p*n₂*n₂ * (∂uex(x_)[1,1]*n¹ + ∂uex(x_)[1,2]*n²)*𝑤
        f2[5:9:27] .+= p*n₂*n₂ * (∂uex(x_)[2,1]*n¹ + ∂uex(x_)[2,2]*n²)*𝑤
        f2[6:9:27] .+= p*n₂*n₂ * (∂uex(x_)[3,1]*n¹ + ∂uex(x_)[3,2]*n²)*𝑤
        f2[7:9:27] .+= 2*p*n₁*n₂ * (∂uex(x_)[1,1]*n¹ + ∂uex(x_)[1,2]*n²)*𝑤
        f2[8:9:27] .+= 2*p*n₁*n₂ * (∂uex(x_)[2,1]*n¹ + ∂uex(x_)[2,2]*n²)*𝑤
        f2[9:9:27] .+= 2*p*n₁*n₂ * (∂uex(x_)[3,1]*n¹ + ∂uex(x_)[3,2]*n²)*𝑤
        f2[1:9:27] .+= p*s₁*n₁ * uex(x_)[1]*Δ
        f2[2:9:27] .+= p*s₁*n₁ * uex(x_)[2]*Δ
        f2[3:9:27] .+= p*s₁*n₁ * uex(x_)[3]*Δ
        f2[4:9:27] .+= p*s₂*n₂ * uex(x_)[1]*Δ
        f2[5:9:27] .+= p*s₂*n₂ * uex(x_)[2]*Δ
        f2[6:9:27] .+= p*s₂*n₂ * uex(x_)[3]*Δ
        f2[7:9:27] .+= p*(s₁*n₂ + s₂*n₁) * uex(x_)[1]*Δ
        f2[8:9:27] .+= p*(s₁*n₂ + s₂*n₁) * uex(x_)[2]*Δ
        f2[9:9:27] .+= p*(s₁*n₂ + s₂*n₁) * uex(x_)[3]*Δ
    end
end
# println(temp)
for a in elements["Ωₚ"][1:1]
# for a in elements["Ωₚ"]
    𝓒 = a.𝓒
    𝓖 = a.𝓖
    for ξ in 𝓖
        x_ = Vec{3}((ξ.x,ξ.y,0.0))
        𝑤 = ξ.𝑤
        Γ¹₁₁ = ξ.Γ¹₁₁
        Γ¹₁₂ = ξ.Γ¹₁₂
        Γ¹₂₂ = ξ.Γ¹₂₂
        Γ²₁₁ = ξ.Γ²₁₁
        Γ²₁₂ = ξ.Γ²₁₂
        Γ²₂₂ = ξ.Γ²₂₂
        Γ¹₁₁₁ = ξ.Γ¹₁₁₁
        Γ¹₁₁₂ = ξ.Γ¹₁₁₂
        Γ¹₁₂₁ = ξ.Γ¹₁₂₁
        Γ¹₁₂₂ = ξ.Γ¹₁₂₂
        Γ¹₂₂₁ = ξ.Γ¹₂₂₁
        Γ¹₂₂₂ = ξ.Γ¹₂₂₂
        Γ²₁₁₁ = ξ.Γ²₁₁₁
        Γ²₁₁₂ = ξ.Γ²₁₁₂
        Γ²₁₂₁ = ξ.Γ²₁₂₁
        Γ²₁₂₂ = ξ.Γ²₁₂₂
        Γ²₂₂₁ = ξ.Γ²₂₂₁
        Γ²₂₂₂ = ξ.Γ²₂₂₂

        Γᵞᵧ₁₁ = Γ¹₁₁₁+Γ²₁₂₁
        Γᵞᵧ₂₂ = Γ¹₁₂₂+Γ²₂₂₂
        Γᵞᵧ₁₂ = Γ¹₁₁₂+Γ²₁₂₂
        Γᵞᵧ₂₁ = Γ¹₁₂₁+Γ²₂₂₁
        Γᵞ₁₁ᵧ = Γ¹₁₁₁+Γ²₁₁₂
        Γᵞ₂₂ᵧ = Γ¹₂₂₁+Γ²₂₂₂
        Γᵞ₁₂ᵧ = Γ¹₁₂₁+Γ²₁₂₂

        Γᵞᵧ₁ = Γ¹₁₁+Γ²₁₂
        Γᵞᵧ₂ = Γ¹₁₂+Γ²₂₂
        Γᵞᵧ₁Γᵝᵦ₁ = Γᵞᵧ₁^2
        Γᵞᵧ₂Γᵝᵦ₂ = Γᵞᵧ₂^2
        Γᵞᵧ₁Γᵝᵦ₂ = Γᵞᵧ₁*Γᵞᵧ₂
        Γᵞᵧ₂Γᵝᵦ₁ = Γᵞᵧ₁*Γᵞᵧ₂
        ΓᵞᵧᵦΓᵝ₁₁ = Γᵞᵧ₁*Γ¹₁₁ + Γᵞᵧ₂*Γ²₁₁
        ΓᵞᵧᵦΓᵝ₂₂ = Γᵞᵧ₁*Γ¹₂₂ + Γᵞᵧ₂*Γ²₂₂
        ΓᵞᵧᵦΓᵝ₁₂ = Γᵞᵧ₁*Γ¹₁₂ + Γᵞᵧ₂*Γ²₁₂


        p = [1,ξ.x,ξ.y]
        p₁ = [0.0,1.0,0.0]
        p₂ = [0.0,0.0,1.0]
        # f2[1:2:6] .+= (p₁ + Γᵞᵧ₁.*p)*𝑤
        # f2[2:2:6] .+= (p₂ + Γᵞᵧ₂.*p)*𝑤
        # f2[1:2:6] .+= p₁*𝑤
        # f2[2:2:6] .+= p₂*𝑤
        # println(Γᵞ₁₁ᵧ*p)
        # println(Γ¹₁₁*p₁+Γ²₁₁*p₂)
        f3[1:9:27] .+= (Γᵞ₁₁ᵧ*p + (Γ¹₁₁*p₁+Γ²₁₁*p₂) + ΓᵞᵧᵦΓᵝ₁₁*p + Γᵞᵧ₁₁*p + 2*Γᵞᵧ₁*p₁ + Γᵞᵧ₁Γᵝᵦ₂*p)*uex(x_)[1]*𝑤
        f3[2:9:27] .+= (Γᵞ₁₁ᵧ*p + (Γ¹₁₁*p₁+Γ²₁₁*p₂) + ΓᵞᵧᵦΓᵝ₁₁*p + Γᵞᵧ₁₁*p + 2*Γᵞᵧ₁*p₁ + Γᵞᵧ₁Γᵝᵦ₂*p)*uex(x_)[2]*𝑤
        f3[3:9:27] .+= (Γᵞ₁₁ᵧ*p + (Γ¹₁₁*p₁+Γ²₁₁*p₂) + ΓᵞᵧᵦΓᵝ₁₁*p + Γᵞᵧ₁₁*p + 2*Γᵞᵧ₁*p₁ + Γᵞᵧ₁Γᵝᵦ₂*p)*uex(x_)[3]*𝑤
        f3[4:9:27] .+= (Γᵞ₂₂ᵧ*p + (Γ¹₂₂*p₁+Γ²₂₂*p₂) + ΓᵞᵧᵦΓᵝ₂₂*p + Γᵞᵧ₂₂*p + 2*Γᵞᵧ₂*p₂ + Γᵞᵧ₂Γᵝᵦ₂*p)*uex(x_)[1]*𝑤
        f3[5:9:27] .+= (Γᵞ₂₂ᵧ*p + (Γ¹₂₂*p₁+Γ²₂₂*p₂) + ΓᵞᵧᵦΓᵝ₂₂*p + Γᵞᵧ₂₂*p + 2*Γᵞᵧ₂*p₂ + Γᵞᵧ₂Γᵝᵦ₂*p)*uex(x_)[2]*𝑤
        f3[6:9:27] .+= (Γᵞ₂₂ᵧ*p + (Γ¹₂₂*p₁+Γ²₂₂*p₂) + ΓᵞᵧᵦΓᵝ₂₂*p + Γᵞᵧ₂₂*p + 2*Γᵞᵧ₂*p₂ + Γᵞᵧ₂Γᵝᵦ₂*p)*uex(x_)[3]*𝑤
        f3[7:9:27] .+= (2*Γᵞ₁₂ᵧ*p + 2*(Γ¹₁₂*p₁+Γ²₁₂*p₂) + 2*ΓᵞᵧᵦΓᵝ₁₂*p + (Γᵞᵧ₁₂ + Γᵞᵧ₂₁)*p + 2*(Γᵞᵧ₁*p₂ + Γᵞᵧ₂*p₁) + 2*Γᵞᵧ₁Γᵝᵦ₂*p)*uex(x_)[1]*𝑤
        f3[8:9:27] .+= (2*Γᵞ₁₂ᵧ*p + 2*(Γ¹₁₂*p₁+Γ²₁₂*p₂) + 2*ΓᵞᵧᵦΓᵝ₁₂*p + (Γᵞᵧ₁₂ + Γᵞᵧ₂₁)*p + 2*(Γᵞᵧ₁*p₂ + Γᵞᵧ₂*p₁) + 2*Γᵞᵧ₁Γᵝᵦ₂*p)*uex(x_)[2]*𝑤
        f3[9:9:27] .+= (2*Γᵞ₁₂ᵧ*p + 2*(Γ¹₁₂*p₁+Γ²₁₂*p₂) + 2*ΓᵞᵧᵦΓᵝ₁₂*p + (Γᵞᵧ₁₂ + Γᵞᵧ₂₁)*p + 2*(Γᵞᵧ₁*p₂ + Γᵞᵧ₂*p₁) + 2*Γᵞᵧ₁Γᵝᵦ₂*p)*uex(x_)[3]*𝑤
        # f3[1:3:9] .+= ((Γ²₁₁*p₁+Γ¹₁₁*p₂))*uex(x_)[1]*𝑤
        # f3[1:3:9] .+= (Γᵞ₁₁ᵧ*p)*uex(x_)[1]*𝑤
    end
end
# println(norm(f1[1:27]-f2))
println(f2)
println(f3)
# println(norm(f1[1:27]-f3))
println(norm(f2+f3))