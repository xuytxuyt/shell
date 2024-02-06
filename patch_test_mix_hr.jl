
using ApproxOperator, Tensors, JLD

import BenchmarkExample: BenchmarkExample
include("import_prescrible_ops.jl")
include("import_patch_test.jl")
elements, nodes = import_roof_mix("msh/patchtest.msh")

𝑅 = BenchmarkExample.PatchTestThinShell.𝑅
𝐿 = BenchmarkExample.PatchTestThinShell.𝐿
E = BenchmarkExample.PatchTestThinShell.𝐸
ν = BenchmarkExample.PatchTestThinShell.𝜈
h = 1.0

cs = BenchmarkExample.cylindricalCoordinate(𝑅)

nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])*3)
s = 3.5*1.0*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
set∇²𝝭!(elements["Ωₚ"])
set∇𝝭!(elements["Γₚ"])
set∇𝝭!(elements["Ωₘ"])
set∇𝝭!(elements["Γₘ"])
set∇𝝭!(elements["Γ¹"])
set∇𝝭!(elements["Γ²"])
set∇𝝭!(elements["Γ³"])
set∇𝝭!(elements["Γ⁴"])
set∇𝝭!(elements["Γ¹ₚ"])
set∇𝝭!(elements["Γ²ₚ"])
set∇𝝭!(elements["Γ³ₚ"])
set∇𝝭!(elements["Γ⁴ₚ"])

n = 1
u(x) = Vec{3}(((x[1]+2*x[2])^n,(3*x[1]+4*x[2])^n,(5*x[1]+6*x[2])^n))
vs = BenchmarkExample.PatchTestThinShell.variables(cs,u)

eval(prescribleForMix)
eval(prescribleBoundary)

eval(opsMix)

f = zeros(3*nₚ)

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
opsh[1](elements["Γ¹ₚ"],elements["Γ¹"],kᴺᵛ,fᴺ)
opsh[1](elements["Γ²ₚ"],elements["Γ²"],kᴺᵛ,fᴺ)
opsh[1](elements["Γ³ₚ"],elements["Γ³"],kᴺᵛ,fᴺ)
opsh[1](elements["Γ⁴ₚ"],elements["Γ⁴"],kᴺᵛ,fᴺ)
opsh[2](elements["Γ¹ₚ"],elements["Γ¹"],kᴹᵛ,fᴹ)
opsh[2](elements["Γ²ₚ"],elements["Γ²"],kᴹᵛ,fᴹ)
opsh[2](elements["Γ³ₚ"],elements["Γ³"],kᴹᵛ,fᴹ)
opsh[2](elements["Γ⁴ₚ"],elements["Γ⁴"],kᴹᵛ,fᴹ)
opsh[3](elements["Γ¹ₚ"],elements["Γ¹"],kᴹᵛ,fᴹ)
opsh[3](elements["Γ²ₚ"],elements["Γ²"],kᴹᵛ,fᴹ)
opsh[3](elements["Γ³ₚ"],elements["Γ³"],kᴹᵛ,fᴹ)
opsh[3](elements["Γ⁴ₚ"],elements["Γ⁴"],kᴹᵛ,fᴹ)

kᵋᵛ = kᴺᵋ\kᴺᵛ
kᵏᵛ = kᴹᵏ\kᴹᵛ
d = (kᵋᵛ'*kᵋᵋ*kᵋᵛ + kᵏᵛ'*kᵏᵏ*kᵏᵛ)\(-f + kᵋᵛ'*kᵋᵋ*(kᴺᵋ\fᴺ) + kᵏᵛ'*kᵏᵏ*(kᴹᵏ\fᴹ))

d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]