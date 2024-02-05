using ApproxOperator, JLD

import BenchmarkExample: BenchmarkExample
include("import_prescrible_ops.jl")
include("import_Scordelis_Lo_roof.jl")
ndiv = 16
elements, nodes = import_roof_mix("msh/scordelislo_"*string(ndiv)*".msh");

𝑅 = BenchmarkExample.ScordelisLoRoof.𝑅
𝐿 = BenchmarkExample.ScordelisLoRoof.𝐿
b₃ = BenchmarkExample.ScordelisLoRoof.𝑞
E = BenchmarkExample.ScordelisLoRoof.𝐸
ν = BenchmarkExample.ScordelisLoRoof.𝜈
h = BenchmarkExample.ScordelisLoRoof.ℎ
cs = BenchmarkExample.cylindricalCoordinate(𝑅)
nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])*3)
s = 3.5*𝐿/2/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ωₘ"])
set∇²𝝭!(elements["Ωₚ"])
set∇𝝭!(elements["Γₘ"])
set∇𝝭!(elements["Γₚ"])
set∇̂³𝝭!(elements["Γᵇ"])
set∇𝝭!(elements["Γʳ"])
set∇̂³𝝭!(elements["Γᵗ"])
set∇̂³𝝭!(elements["Γˡ"])
set𝝭!(elements["𝐴"])

eval(prescribleForMix)
eval(prescribleBoundary)

eval(opsMix)

opForce = Operator{:∫vᵢbᵢdΩ}()

f = zeros(3*nₚ)
opForce(elements["Ωₘ"],f)

kᵋᵋ = zeros(6*nᵥ,6*nᵥ)
kᴺᵋ = zeros(6*nᵥ,6*nᵥ)
kᴺᵛ = zeros(6*nᵥ,3*nₚ)
kᵏᵏ = zeros(9*nᵥ,9*nᵥ)
kᴹᵏ = zeros(9*nᵥ,9*nᵥ)
kᴹᵛ = zeros(9*nᵥ,3*nₚ)
fᴺ = zeros(6*nᵥ)
fᴹ = zeros(3*nᵥ)

ops[1](elements["Ω"],kᵋᵋ)
ops[2](elements["Ω"],kᴺᵋ)
ops[3](elements["Γₚ"],elements["Γₘ"],kᴺᵛ)
ops[4](elements["Ωₚ"],elements["Ωₘ"],kᴺᵛ)

ops[5](elements["Ω"],kᵏᵏ)
ops[6](elements["Ω"],kᴹᵏ)
ops[7](elements["Γₚ"],elements["Γₘ"],kᴹᵛ)
ops[8](elements["Γₚ"],elements["Γₘ"],kᴹᵛ)
ops[9](elements["Γₚ"],elements["Γₘ"],kᴹᵛ)
ops[10](elements["Ωₚ"],elements["Ωₘ"],kᴹᵛ)

eval(opsNitsche)
kᵛ = zeros(3*nₚ,3*nₚ)
fᵛ = zeros(3*nₚ)
opsv[1](elements["Γᵇ"],kᵛ,fᵛ)
opsv[1](elements["Γᵗ"],kᵛ,fᵛ)
opsv[1](elements["Γˡ"],kᵛ,fᵛ)
opsv[2](elements["Γᵇ"],kᵛ,fᵛ)
opsv[2](elements["Γᵗ"],kᵛ,fᵛ)
opsv[2](elements["Γˡ"],kᵛ,fᵛ)
opsv[3](elements["Γᵗ"],kᵛ,fᵛ)
opsv[3](elements["Γˡ"],kᵛ,fᵛ)

αᵥ = 1e5
αᵣ = 1e3
eval(opsPenalty)
kᵅ = zeros(3*nₚ,3*nₚ)
fᵅ = zeros(3*nₚ)
opsα[1](elements["Γᵇ"],kᵅ,fᵅ)
opsα[1](elements["Γᵗ"],kᵅ,fᵅ)
opsα[1](elements["Γˡ"],kᵅ,fᵅ)
opsα[2](elements["Γᵗ"],kᵅ,fᵅ)
opsα[2](elements["Γˡ"],kᵅ,fᵅ)

d = (kᵛ+kᵅ + (kᴺᵋ\kᴺᵛ)'*kᵋᵋ*(kᴺᵋ\kᴺᵛ) + (kᴹᵏ\kᴹᵛ)'*kᵏᵏ*(kᴹᵏ\kᴹᵛ))\(f+fᵛ+fᵅ)

d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

op𝐴 = Operator{:ScordelisLoRoof_𝐴}()
push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
w = op𝐴(elements["𝐴"])

println(w)
@save compress=true "jld/scordelislo_mix_nitsche_"*string(ndiv)*".jld" d₁ d₂ d₃