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

eval(prescribleForMix)
eval(prescribleBoundary)

ops = [
    Operator{:∫NC⁻¹NdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫MC⁻¹MdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫𝒏𝑵𝒗dΓ}(),
    Operator{:∫∇𝑵𝒗dΩ}(),
    Operator{:∫∇𝑴𝒏𝒂₃𝒗dΓ}(),
    Operator{:∫𝑴ₙₙ𝜽ₙdΓ}(),
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

# ops[1](elements["Ωₚ"][1:1],kᴺᴺ)
ops[1](elements["Ωₚ"],kᴺᴺ)
ops[2](elements["Ωₚ"],kᴹᴹ)
ops[3](elements["Γₚ"],elements["Γ"],kᴺᵛ)
ops[4](elements["Ωₚ"],elements["Ω"],kᴺᵛ)
ops[5](elements["Γₚ"],elements["Γ"],kᴹᵛ)
ops[6](elements["Γₚ"],elements["Γ"],kᴹᵛ)
ops[7](elements["Γₚ"],elements["Γ"],kᴹᵛ)
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

# k = [kᴹᴹ kᴹᵛ;kᴹᵛ' kᵛᵛ]
# f = [fᴹ;fᵛ]
# d = k\f
# d₁ = d[9*nᵥ+1:3:end]
# d₂ = d[9*nᵥ+2:3:end]
# d₃ = d[9*nᵥ+3:3:end]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
w = ops[12](elements["𝐴"])

println(w)
# @save compress=true "jld/scordelislo_mix_"*string(ndiv)*".jld" d₁ d₂ d₃