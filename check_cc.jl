
using ApproxOperator, Tensors, BenchmarkExample, LinearAlgebra

include("import_Scordelis_Lo_roof.jl")
ndiv = 16
elements, nodes = import_roof_mix("msh/scordelislo_"*string(ndiv)*".msh");

𝑅 = BenchmarkExample.ScordelisLoRoof.𝑅
𝐿 = BenchmarkExample.ScordelisLoRoof.𝐿
b₃ = BenchmarkExample.ScordelisLoRoof.𝑞
E = BenchmarkExample.ScordelisLoRoof.𝐸
ν = BenchmarkExample.ScordelisLoRoof.𝜈
h = BenchmarkExample.ScordelisLoRoof.ℎ
cs = BenchmarkExample.cylindricalCoordinate(BenchmarkExample.ScordelisLoRoof.𝑅)

nₚ = length(nodes)
nₑ = length(elements["Ω"])
nᵥ = nₑ*3

s = 3.5*𝐿/2/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set∇𝝭!(elements["Ωₘ"])
set∇²𝝭!(elements["Ωₚ"])
set∇𝝭!(elements["Γₘ"])
set∇𝝭!(elements["Γₚ"])

eval(prescribleForMix)
eval(prescribleBoundary)

ops = [
    Operator{:∫δεCεdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫δNεdΩ}(),
    Operator{:∫𝒏𝑵𝒗dΓ}(),
    Operator{:∫∇𝑵𝒗dΩ}(),
    Operator{:∫δκCκdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫δMκdΩ}(),
    Operator{:∫∇𝑴𝒏𝒂₃𝒗dΓ}(),
    Operator{:∫𝑴ₙₙ𝜽ₙdΓ}(),
    Operator{:ΔMₙₛ𝒂₃𝒗}(),
    Operator{:∫∇𝑴∇𝒂₃𝒗dΩ}(),
]

kᵋᵋ = zeros(6*nᵥ,6*nᵥ)
kᴺᵋ = zeros(6*nᵥ,6*nᵥ)
kᴺᵛ = zeros(6*nᵥ,3*nₚ)
kᵏᵏ = zeros(9*nᵥ,9*nᵥ)
kᴹᵏ = zeros(9*nᵥ,9*nᵥ)
kᴹᵛ = zeros(9*nᵥ,3*nₚ)

ops[1](elements["Ωₚ"],kᵋᵋ)
ops[2](elements["Ωₚ"],kᴺᵋ)
ops[3](elements["Γₚ"],elements["Γₘ"],kᴺᵛ)
ops[4](elements["Ωₚ"],elements["Ωₘ"],kᴺᵛ)

ops[5](elements["Ωₚ"],kᵏᵏ)
ops[6](elements["Ωₚ"],kᴹᵏ)
ops[7](elements["Γₚ"],elements["Γₘ"],kᴹᵛ)
ops[8](elements["Γₚ"],elements["Γₘ"],kᴹᵛ)
ops[9](elements["Γₚ"],elements["Γₘ"],kᴹᵛ)
ops[10](elements["Ωₚ"],elements["Ωₘ"],kᴹᵛ)

kᵇ = - kᴺᵋ\kᴺᵛ
kᵐ = - kᴹᵏ\kᴹᵛ

n = 1
uex(x) = Vec{3}(((x[1]+2*x[2])^n,(3*x[1]+4*x[2])^n,(5*x[1]+6*x[2])^n))
∂uex(x) = gradient(uex,x)
∂²uex(x) = gradient(∂uex,x)
dᵛ = zeros(3*nₚ)
for (I,node) in enumerate(nodes)
    x = Vec{3}((node.x,node.y,node.z))
    u = uex(x)
    dᵛ[3*I-2] = u[1]
    dᵛ[3*I-1] = u[2]
    dᵛ[3*I]   = u[3]
end

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
        # println(∂₁u[1])
        # println(∂₁₁uʰ[1])
        # global err += abs(∂₁u[1]-∂₁uʰ[1])
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
# # println(f)