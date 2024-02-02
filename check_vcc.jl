
using ApproxOperator, Tensors

include("import_Scordelis_Lo_roof.jl")
ndiv = 16
elements, nodes = import_roof_mix("msh/scordelislo_"*string(ndiv)*".msh",ndiv-1);

nₚ = length(nodes)
nₑ = length(elements["Ω"])
nᵥ = nₑ*3

s = 3.5*𝐿/2/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set∇²𝝭!(elements["Ω"])
set∇²𝝭!(elements["Ωₚ"])
set∇𝝭!(elements["Γ"])
set∇𝝭!(elements["Γₚ"])

eval(prescribleForMix)
eval(prescribleBoundary)

kᴺᴺ = zeros(3*nᵥ,3*nᵥ)
kᴺᵛ = zeros(3*nᵥ,3*nₚ)

ops = [
    Operator{:∫NC⁻¹NdΩ}(:E=>1.0,:ν=>0.0,:h=>0.0),
    Operator{:∫MC⁻¹MdΩ}(:E=>1.0,:ν=>0.0,:h=>0.0),
    Operator{:∫𝒏𝑵𝒗dΓ}(),
    Operator{:∫∇𝑵𝒗dΩ}(),
    Operator{:∫∇𝑴𝒏𝒂₃𝒗dΓ}(),
    Operator{:∫𝑴ₙₙ𝜽ₙdΓ}(),
    Operator{:ΔMₙₛ𝒂₃𝒗}(),
    Operator{:∫∇𝑴∇𝒂₃𝒗dΩ}()
]

ops[1](elements["Ωₚ"],kᴺᴺ)
ops[3](elements["Γₚ"],elements["Γ"],kᴺᵛ)
# ops[4](elements["Ωₚ"],elements["Ω"],kᴺᵛ)

uex(x) = Vec{3}((1.0,1.0,1.0))
dᵛ = zeros(3*nₚ)
for (I,node) in enumerate(nodes)
    x = Vec{3}((node.x,node.y,node.z))
    u = uex(x)
    dᵛ[3*I-2] = u[1]
    dᵛ[3*I-1] = u[2]
    dᵛ[3*I]   = u[3]
end
fᴺ = kᴺᵛ*dᵛ

fᴺ_ = zeros(3*nᵥ)
for a in elements["Γₚ"]
    𝓒 = a.𝓒
    𝓖 = a.𝓖
    for ξ in 𝓖
        𝑤 = ξ.𝑤
        N = ξ[:𝝭]
        n₁ = ξ.n₁
        n₂ = ξ.n₂
        𝒂₁₍₁₎ = ξ.𝒂₁₍₁₎
        𝒂₁₍₂₎ = ξ.𝒂₁₍₂₎
        𝒂₁₍₃₎ = ξ.𝒂₁₍₃₎
        𝒂₂₍₁₎ = ξ.𝒂₂₍₁₎
        𝒂₂₍₂₎ = ξ.𝒂₂₍₂₎
        𝒂₂₍₃₎ = ξ.𝒂₂₍₃₎
        x = Vec{3}((ξ.x,ξ.y,ξ.z))
        𝒂₁ = Vec{3}((𝒂₁₍₁₎,𝒂₁₍₂₎,𝒂₁₍₃₎))
        𝒂₂ = Vec{3}((𝒂₂₍₁₎,𝒂₂₍₂₎,𝒂₂₍₃₎))
        v₁ = 𝒂₁⋅uex(x)
        v₂ = 𝒂₂⋅uex(x)
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            fᴺ_[3*I-2] -= N[i]*v₁*n₁*𝑤
            fᴺ_[3*I-1] -= N[i]*v₂*n₂*𝑤
            fᴺ_[3*I]   -= N[i]*(v₁*n₂+v₂*n₁)*𝑤
        end
    end
end