using ApproxOperator, JLD

import BenchmarkExample: BenchmarkExample
include("import_Scordelis_Lo_roof.jl")
ndiv = 11
elements, nodes = import_roof_mix("msh/scordelislo_"*string(ndiv)*".msh",ndiv-1);

𝑅 = BenchmarkExample.ScordelisLoRoof.𝑅
𝐿 = BenchmarkExample.ScordelisLoRoof.𝐿
b₃ = BenchmarkExample.ScordelisLoRoof.𝑞
E = BenchmarkExample.ScordelisLoRoof.𝐸
ν = BenchmarkExample.ScordelisLoRoof.𝜈
h = BenchmarkExample.ScordelisLoRoof.ℎ
cs = BenchmarkExample.cylindricalCoordinate(𝑅)

nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])/8*3)
s = 2.5*𝐿/2/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set∇𝝭!(elements["Ω"])
set∇²𝝭!(elements["Ωₚ"])
set∇𝝭!(elements["Γ"])
set∇𝝭!(elements["Γₚ"])
set𝝭!(elements["Γᵇ"])
set𝝭!(elements["Γʳ"])
set∇𝝭!(elements["Γᵗ"])
set∇𝝭!(elements["Γˡ"])
set𝝭!(elements["𝐴"])

prescribe!(elements["Ωₚ"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:a²²=>(ξ¹,ξ²,ξ³)->cs.a²²(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:a¹²=>(ξ¹,ξ²,ξ³)->cs.a¹²(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:𝒂¹₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(ξ¹,ξ²)[1])
prescribe!(elements["Ωₚ"],:𝒂¹₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(ξ¹,ξ²)[2])
prescribe!(elements["Ωₚ"],:𝒂¹₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(ξ¹,ξ²)[3])
prescribe!(elements["Ωₚ"],:𝒂²₍₁₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(ξ¹,ξ²)[1])
prescribe!(elements["Ωₚ"],:𝒂²₍₂₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(ξ¹,ξ²)[2])
prescribe!(elements["Ωₚ"],:𝒂²₍₃₎=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(ξ¹,ξ²)[3])
prescribe!(elements["Ωₚ"],:Γ¹₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ¹₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ¹₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₂₂(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ²₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₁(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ²₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₂(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ²₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ¹₁₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁₁(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ¹₁₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁₂(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ¹₁₂₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂₁(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ¹₁₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂₂(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ¹₂₂₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₂₂₁(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ²₁₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₁₂(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ²₁₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₂₂(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ²₂₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₁₁(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ²₂₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₁₂(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ²₂₂₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂₁(ξ¹,ξ²))
prescribe!(elements["Ωₚ"],:Γ²₂₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂₂(ξ¹,ξ²))

prescribe!(elements["Γₚ"],:Γ¹₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁(ξ¹,ξ²))
prescribe!(elements["Γₚ"],:Γ¹₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂(ξ¹,ξ²))
prescribe!(elements["Γₚ"],:Γ¹₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₂₂(ξ¹,ξ²))
prescribe!(elements["Γₚ"],:Γ²₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₁(ξ¹,ξ²))
prescribe!(elements["Γₚ"],:Γ²₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₂(ξ¹,ξ²))
prescribe!(elements["Γₚ"],:Γ²₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂(ξ¹,ξ²))
prescribe!(elements["Γₚ"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(ξ¹,ξ²))
prescribe!(elements["Γₚ"],:a²²=>(ξ¹,ξ²,ξ³)->cs.a²²(ξ¹,ξ²))
prescribe!(elements["Γₚ"],:a¹²=>(ξ¹,ξ²,ξ³)->cs.a¹²(ξ¹,ξ²))

prescribe!(elements["Γ"],:n₁₁=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γ"],:n₂₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γ"],:n₃₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γ"],:θ=>(ξ¹,ξ²,ξ³)->0.0)

prescribe!(elements["Γᵇ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:n₁₁=>(ξ¹,ξ²,ξ³)->1.0)
prescribe!(elements["Γᵇ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:n₂₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:n₃₃=>(ξ¹,ξ²,ξ³)->1.0)
prescribe!(elements["Γᵗ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:n₁₁=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:n₂₂=>(ξ¹,ξ²,ξ³)->1.0)
prescribe!(elements["Γᵗ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:n₃₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:a₃₁=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(ξ¹,ξ²)[1])
prescribe!(elements["Γᵗ"],:a₃₂=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(ξ¹,ξ²)[2])
prescribe!(elements["Γᵗ"],:a₃₃=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(ξ¹,ξ²)[3])
prescribe!(elements["Γᵗ"],:θ=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:n₁₁=>(ξ¹,ξ²,ξ³)->1.0)
prescribe!(elements["Γˡ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:n₂₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:n₃₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:θ=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:a₃₁=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(ξ¹,ξ²)[1])
prescribe!(elements["Γˡ"],:a₃₂=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(ξ¹,ξ²)[2])
prescribe!(elements["Γˡ"],:a₃₃=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(ξ¹,ξ²)[3])
prescribe!(elements["Ω"],:b₁=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Ω"],:b₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Ω"],:b₃=>(ξ¹,ξ²,ξ³)->b₃)

ops = [
    Operator{:∫NC⁻¹NdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫MC⁻¹MdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫𝒏𝑵𝒗dΓ}(),
    Operator{:∫∇𝑵𝒗dΩ}(),
    Operator{:∫∇𝑴𝒏𝒂₃𝒗dΓ}(),
    Operator{:∫MₙₙθₙdΓ}(),
    Operator{:ΔMₙₛ𝒂₃𝒗}(),
    Operator{:∫∇𝑴∇𝒂₃𝒗dΩ}(),
    Operator{:∫vᵢbᵢdΩ}(),
    Operator{:∫vᵢgᵢdΓ}(:α=>1e9*E),
    Operator{:∫δθθdΓ}(:α=>1e7*E),
    Operator{:ScordelisLoRoof_𝐴}()
]
kᴺᴺ = zeros(6*nᵥ,6*nᵥ)
kᴹᴹ = zeros(9*nᵥ,9*nᵥ)
kᴺᴹ = zeros(6*nᵥ,9*nᵥ)
kᴺᵛ = zeros(6*nᵥ,3*nₚ)
kᴹᵛ = zeros(9*nᵥ,3*nₚ)
kᵛᵛ = zeros(3*nₚ,3*nₚ)
fᴺ = zeros(6*nᵥ)
fᴹ = zeros(9*nᵥ)
fᵛ = zeros(3*nₚ)

ops[1](elements["Ωₚ"],kᴺᴺ)
ops[2](elements["Ωₚ"],kᴹᴹ)
ops[3](elements["Γₚ"],elements["Γ"],kᴺᵛ,fᴺ)
ops[4](elements["Ωₚ"],elements["Ω"],kᴺᵛ)
ops[5](elements["Γₚ"],elements["Γ"],kᴹᵛ,fᴹ)
ops[6](elements["Γₚ"],elements["Γ"],kᴹᵛ,fᴹ)
ops[7](elements["Γₚ"],elements["Γ"],kᴹᵛ,fᴹ)
ops[8](elements["Ωₚ"],elements["Ω"],kᴹᵛ)

ops[9](elements["Ω"],fᵛ)
ops[10](elements["Γᵇ"],kᵛᵛ,fᵛ)
ops[10](elements["Γᵗ"],kᵛᵛ,fᵛ)
ops[10](elements["Γˡ"],kᵛᵛ,fᵛ)
ops[11](elements["Γᵗ"],kᵛᵛ,fᵛ)
ops[11](elements["Γˡ"],kᵛᵛ,fᵛ)

k = [kᴺᴺ kᴺᴹ kᴺᵛ;kᴺᴹ' kᴹᴹ kᴹᵛ;kᴺᵛ' kᴹᵛ' kᵛᵛ]
f = [fᴺ;fᴹ;fᵛ]
d = k\f
d₁ = d[15*nᵥ+1:3:end]
d₂ = d[15*nᵥ+2:3:end]
d₃ = d[15*nᵥ+3:3:end]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
w = ops[12](elements["𝐴"])

# @save compress=true "jld/scordelislo_gauss_"*string(ndiv)*".jld" d₁ d₂ d₃