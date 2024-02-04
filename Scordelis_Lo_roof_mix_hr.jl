using ApproxOperator, JLD

import BenchmarkExample: BenchmarkExample
include("import_Scordelis_Lo_roof.jl")
ndiv = 16
elements, nodes = import_roof_mix("msh/scordelislo_"*string(ndiv)*".msh",ndiv-1);

𝑅 = BenchmarkExample.ScordelisLoRoof.𝑅
𝐿 = BenchmarkExample.ScordelisLoRoof.𝐿
b₃ = BenchmarkExample.ScordelisLoRoof.𝑞
E = BenchmarkExample.ScordelisLoRoof.𝐸
ν = BenchmarkExample.ScordelisLoRoof.𝜈
h = BenchmarkExample.ScordelisLoRoof.ℎ
cs = BenchmarkExample.cylindricalCoordinate(𝑅)
nₚ = length(nodes)
# nᵥ = Int(length(elements["Ω"])/2*3)
nᵥ = Int(length(elements["Ω"])*3)
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
    Operator{:ScordelisLoRoof_𝐴}()
]
opForce = Operator{:∫vᵢbᵢdΩ}()

opΩ = [
    Operator{:∫εᵢⱼNᵢⱼdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫κᵢⱼMᵢⱼdΩ}(:E=>E,:ν=>ν,:h=>h),
]

opPenalty = [
    Operator{:∫vᵢgᵢdΓ}(:α=>1e9*E),
    Operator{:∫δθθdΓ}(:α=>1e7*E),
]

kᵋᵋ = zeros(6*nᵥ,6*nᵥ)
kᴺᵋ = zeros(6*nᵥ,6*nᵥ)
kᴺᵛ = zeros(6*nᵥ,3*nₚ)
kᵏᵏ = zeros(9*nᵥ,9*nᵥ)
kᴹᵏ = zeros(9*nᵥ,9*nᵥ)
kᴹᵛ = zeros(9*nᵥ,3*nₚ)
fᴺ = zeros(6*nᵥ)
fᴹ = zeros(3*nᵥ)

ops[1](elements["Ωₚ"],kᵋᵋ)
ops[2](elements["Ωₚ"],kᴺᵋ)
ops[3](elements["Γₚ"],elements["Γ"],kᴺᵛ)
ops[4](elements["Ωₚ"],elements["Ω"],kᴺᵛ)

ops[5](elements["Ωₚ"],kᵏᵏ)
ops[6](elements["Ωₚ"],kᴹᵏ)
ops[7](elements["Γₚ"],elements["Γ"],kᴹᵛ)
ops[8](elements["Γₚ"],elements["Γ"],kᴹᵛ)
ops[9](elements["Γₚ"],elements["Γ"],kᴹᵛ)
ops[10](elements["Ωₚ"],elements["Ω"],kᴹᵛ)

k = zeros(3*nₚ,3*nₚ)
f = zeros(3*nₚ)
kᵅ = zeros(3*nₚ,3*nₚ)
fᵅ = zeros(3*nₚ)
kᵝ = zeros(3*nₚ,3*nₚ)
fᵝ = zeros(3*nₚ)
opForce(elements["Ω"],f)
opPenalty[1](elements["Γᵇ"],kᵅ,fᵅ)
opPenalty[1](elements["Γᵗ"],kᵅ,fᵅ)
opPenalty[1](elements["Γˡ"],kᵅ,fᵅ)
opPenalty[2](elements["Γᵗ"],kᵝ,fᵝ)
opPenalty[2](elements["Γˡ"],kᵝ,fᵝ)

# opΩ[1](elements["Ω"],k)
# opΩ[2](elements["Ω"],k)

# d = [k+kᵅ+kᵝ zeros(3*nₚ,6*nᵥ) kᴺᵛ';zeros(6*nᵥ,3*nₚ) kᵋᵋ kᴺᵋ';kᴺᵛ kᴺᵋ zeros(6*nᵥ,6*nᵥ)]\[f+fᵅ+fᵝ;zeros(12*nᵥ)]
# d = (k+kᵅ+kᵝ + (kᴺᵋ\kᴺᵛ)'*kᵋᵋ*(kᴺᵋ\kᴺᵛ))\(f+fᵅ+fᵝ)
d = (kᵅ+kᵝ + (kᴺᵋ\kᴺᵛ)'*kᵋᵋ*(kᴺᵋ\kᴺᵛ) + (kᴹᵏ\kᴹᵛ)'*kᵏᵏ*(kᴹᵏ\kᴹᵛ))\(f+fᵅ+fᵝ)
# d = (k+kᵅ+kᵝ)\(f+fᵅ+fᵝ)
# d = [k+kᵅ+kᵝ kᴹᵛ';kᴹᵛ kᴹᴹ]\[f+fᵅ+fᵝ;fᴹ]
# d = [k+kᵅ+kᵝ kᴺᵛ';kᴺᵛ kᴺᴺ]\[-f+fᵅ+fᵝ;fᴺ]
# d = [zeros(3*nₚ,3*nₚ) kᴺᵛ' kᴹᵛ';kᴺᵛ kᴺᴺ zeros(3*nᵥ,3*nᵥ);kᴹᵛ zeros(3*nᵥ,3*nᵥ) kᴹᴹ]\[f+fᵅ+fᵝ;fᴺ;fᴹ]
# d = (kᴺᵛ'*(kᴺᴺ\kᴺᵛ) + kᴹᵛ'*(kᴹᴹ\kᴹᵛ) + kᵅ + kᵝ)\(-f+fᵅ+fᵝ)
# d = (kᴺᵛ'*(kᴺᴺ\kᴺᵛ) + k + kᵅ + kᵝ)\(f+fᵅ+fᵝ)
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
w = ops[11](elements["𝐴"])

println(w)
@save compress=true "jld/scordelislo_mix_"*string(ndiv)*".jld" d₁ d₂ d₃