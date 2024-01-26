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
# h = 1.0
nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])/2*3)
s = 3.5*𝐿/2/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set∇²𝝭!(elements["Ω"])
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
opΩ = [
    Operator{:∫εᵢⱼNᵢⱼdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫κᵢⱼMᵢⱼdΩ}(:E=>E,:ν=>ν,:h=>h),
]
kᴺᴺ = zeros(3*nᵥ,3*nᵥ)
kᴹᴹ = zeros(3*nᵥ,3*nᵥ)
kᴺᴹ = zeros(3*nᵥ,3*nᵥ)
kᴺᵛ = zeros(3*nᵥ,3*nₚ)
kᴹᵛ = zeros(3*nᵥ,3*nₚ)
fᴺ = zeros(3*nᵥ)
fᴹ = zeros(3*nᵥ)

ops[1](elements["Ωₚ"],kᴺᴺ)
ops[2](elements["Ωₚ"],kᴹᴹ)
ops[3](elements["Γₚ"],elements["Γ"],kᴺᵛ)
ops[4](elements["Ωₚ"],elements["Ω"],kᴺᵛ)
ops[5](elements["Γₚ"],elements["Γ"],kᴹᵛ)
ops[6](elements["Γₚ"],elements["Γ"],kᴹᵛ)
ops[7](elements["Γₚ"],elements["Γ"],kᴹᵛ)
ops[8](elements["Ωₚ"],elements["Ω"],kᴹᵛ)

k = zeros(3*nₚ,3*nₚ)
f = zeros(3*nₚ)
kᵅ = zeros(3*nₚ,3*nₚ)
fᵅ = zeros(3*nₚ)
kᵝ = zeros(3*nₚ,3*nₚ)
fᵝ = zeros(3*nₚ)
ops[9](elements["Ω"],f)
ops[10](elements["Γᵇ"],kᵅ,fᵅ)
ops[10](elements["Γᵗ"],kᵅ,fᵅ)
ops[10](elements["Γˡ"],kᵅ,fᵅ)
ops[11](elements["Γᵗ"],kᵝ,fᵝ)
ops[11](elements["Γˡ"],kᵝ,fᵝ)

# opΩ[1](elements["Ω"],k)
opΩ[2](elements["Ω"],k)

# d = (k+kᵅ+kᵝ)\(f+fᵅ+fᵝ)
# d = [k+kᵅ+kᵝ kᴹᵛ';kᴹᵛ kᴹᴹ]\[f+fᵅ+fᵝ;fᴹ]
# d = [k+kᵅ+kᵝ kᴺᵛ';kᴺᵛ kᴺᴺ]\[-f+fᵅ+fᵝ;fᴺ]
# d = [zeros(3*nₚ,3*nₚ) kᴺᵛ' kᴹᵛ';kᴺᵛ kᴺᴺ zeros(3*nᵥ,3*nᵥ);kᴹᵛ zeros(3*nᵥ,3*nᵥ) kᴹᴹ]\[f+fᵅ+fᵝ;fᴺ;fᴹ]
# d = (kᴺᵛ'*(kᴺᴺ\kᴺᵛ) + kᴹᵛ'*(kᴹᴹ\kᴹᵛ) + kᵅ + kᵝ)\(f+fᵅ+fᵝ)
d = (kᴺᵛ'*(kᴺᴺ\kᴺᵛ) + k + kᵅ + kᵝ)\(f+fᵅ+fᵝ)
# d = (k + kᵅ + kᵝ)\(f+fᵅ+fᵝ)
d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

# k = [kᴹᴹ kᴹᵛ;kᴹᵛ' kᵛᵛ]
# f = [fᴹ;fᵛ]
# d = k\f
# d₁ = d[9*nᵥ+1:3:end]
# d₂ = d[9*nᵥ+2:3:end]
# d₃ = d[9*nᵥ+3:3:end]

# d = (kᴺᵛ'*(kᴺᴺ\kᴺᵛ) + kᴹᵛ'*(kᴹᴹ\kᴹᵛ) + kᵛᵛ)\fᵛ
# d = (kᴹᵛ'*(kᴹᴹ\kᴹᵛ) + kᵛᵛ)\fᵛ

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
w = ops[12](elements["𝐴"])

println(w)
@save compress=true "jld/scordelislo_mix_"*string(ndiv)*".jld" d₁ d₂ d₃