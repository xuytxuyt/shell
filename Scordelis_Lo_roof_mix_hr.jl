using ApproxOperator, JLD

import BenchmarkExample: BenchmarkExample
include("import_prescrible_ops.jl")
include("import_Scordelis_Lo_roof.jl")
ndiv = 28
elements, nodes = import_roof_mix("msh/scordelislo_"*string(ndiv)*".msh");

nₘ = 21
𝑅 = BenchmarkExample.ScordelisLoRoof.𝑅
𝐿 = BenchmarkExample.ScordelisLoRoof.𝐿
b₃ = BenchmarkExample.ScordelisLoRoof.𝑞
E = BenchmarkExample.ScordelisLoRoof.𝐸
ν = BenchmarkExample.ScordelisLoRoof.𝜈
h = BenchmarkExample.ScordelisLoRoof.ℎ

cs = BenchmarkExample.cylindricalCoordinate(𝑅)
nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])*3)
s = 2.5*𝐿/2/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

eval(prescribeForMix)
eval(prescribeVariables)

set𝝭!(elements["Ω"])
set∇²𝝭!(elements["Ωₚ"])
set∇𝝭!(elements["Γₚ"])
set∇𝝭!(elements["Ωₘ"])
set∇𝝭!(elements["Γₘ"])
set∇𝝭!(elements["Γᵇ"])
set∇𝝭!(elements["Γᵗ"])
set∇𝝭!(elements["Γˡ"])
set𝝭!(elements["𝐴"])

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

eval(opsHR)
fᴺ = zeros(6*nᵥ)
fᴹ = zeros(9*nᵥ)
opsh[1](elements["Γᵇₚ"],elements["Γᵇ"],kᴺᵛ,fᴺ)
opsh[1](elements["Γᵗₚ"],elements["Γᵗ"],kᴺᵛ,fᴺ)
opsh[1](elements["Γˡₚ"],elements["Γˡ"],kᴺᵛ,fᴺ)
opsh[2](elements["Γᵇₚ"],elements["Γᵇ"],kᴹᵛ,fᴹ)
opsh[2](elements["Γᵗₚ"],elements["Γᵗ"],kᴹᵛ,fᴹ)
opsh[2](elements["Γˡₚ"],elements["Γˡ"],kᴹᵛ,fᴹ)
opsh[3](elements["Γᵗₚ"],elements["Γᵗ"],kᴹᵛ,fᴹ)
opsh[3](elements["Γˡₚ"],elements["Γˡ"],kᴹᵛ,fᴹ)

kᵋᵛ = kᴺᵋ\kᴺᵛ
kᵏᵛ = kᴹᵏ\kᴹᵛ
d = (kᵋᵛ'*kᵋᵋ*kᵋᵛ + kᵏᵛ'*kᵏᵏ*kᵏᵛ)\(-f + kᵋᵛ'*kᵋᵋ*(kᴺᵋ\fᴺ) + kᵏᵛ'*kᵏᵏ*(kᴹᵏ\fᴹ))

d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

op𝐴 = Operator{:ScordelisLoRoof_𝐴}()
push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
w = op𝐴(elements["𝐴"])

println(w)
@save compress=true "jld/scordelislo_mix_hr_"*string(ndiv)*".jld" d₁ d₂ d₃