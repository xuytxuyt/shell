
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
# cs = BenchmarkExample.cartesianCoordinate()

nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])*3)
s = 3.5*0.2*ones(nₚ)
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
# u(x) = Vec{3}(((x[1]+2*x[2])^n,(3*x[1]+4*x[2])^n,(5*x[1]+6*x[2])^n))
# u(x) = Vec{3}(((x[1]+2*x[2])^n,(3*x[1]+4*x[2])^n,0.0))
u(x) = Vec{3}((0.0,0.0,(5*x[1]+6*x[2])^n))
# u(x) = Vec{3}((0.0,3*x[1]+4*x[2],0.0))
vs = BenchmarkExample.PatchTestThinShell.variables(cs,u)

eval(prescribeForMix)
eval(prescribeForPenalty)
eval(prescribeVariables)

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

αᵥ = 1e9
αᵣ = 1e7
eval(opsPenalty)
kᵅ = zeros(3*nₚ,3*nₚ)
fᵅ = zeros(3*nₚ)
opsα[1](elements["Γ¹"],kᵅ,fᵅ)
opsα[1](elements["Γ²"],kᵅ,fᵅ)
opsα[1](elements["Γ³"],kᵅ,fᵅ)
opsα[1](elements["Γ⁴"],kᵅ,fᵅ)
opsα[2](elements["Γ¹"],kᵅ,fᵅ)
opsα[2](elements["Γ²"],kᵅ,fᵅ)
opsα[2](elements["Γ³"],kᵅ,fᵅ)
opsα[2](elements["Γ⁴"],kᵅ,fᵅ)

kᵋᵛ = kᴺᵋ\kᴺᵛ
kᵏᵛ = kᴹᵏ\kᴹᵛ
d = (kᵅ + kᵋᵛ'*kᵋᵋ*kᵋᵛ + kᵏᵛ'*kᵏᵏ*kᵏᵛ)\(-f + fᵅ)

d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)

set∇²𝝭!(elements["Ωᵍ"])
opE = Operator{:Hₑ_ThinShell}(:E=>E,:ν=>ν,:h=>h)
Hₑ, L₂ = opE(elements["Ωᵍ"])