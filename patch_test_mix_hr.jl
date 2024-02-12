
using ApproxOperator, Tensors, JLD, LinearAlgebra

import BenchmarkExample: BenchmarkExample
include("import_prescrible_ops.jl")
include("import_patch_test.jl")
elements, nodes = import_roof_mix("msh/patchtest.msh")

# nₘ = 55
nₘ = 21
𝑅 = BenchmarkExample.PatchTestThinShell.𝑅
𝐿 = BenchmarkExample.PatchTestThinShell.𝐿
E = BenchmarkExample.PatchTestThinShell.𝐸
ν = BenchmarkExample.PatchTestThinShell.𝜈
h = BenchmarkExample.PatchTestThinShell.ℎ

cs = BenchmarkExample.cylindricalCoordinate(𝑅)
# cs = BenchmarkExample.cartesianCoordinate()

nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])*3)
s = 2.5*0.1*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

n = 2
u(x) = Vec{3}(((x[1]+2*x[2])^n,(3*x[1]+4*x[2])^n,(5*x[1]+6*x[2])^n))
vs = BenchmarkExample.PatchTestThinShell.variables(cs,u)

eval(prescribeForMix)
eval(prescribeVariables)

set𝝭!(elements["Ω"])
set∇²𝝭!(elements["Ωₚ"])
set∇𝝭!(elements["Γₚ"])
set∇𝝭!(elements["Ωₘ"])
set∇𝝭!(elements["Γₘ"])
set∇𝝭!(elements["Γ¹"])
set∇𝝭!(elements["Γ²"])
set∇𝝭!(elements["Γ³"])
set∇𝝭!(elements["Γ⁴"])

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
opsh[4](elements["Γ¹ₚ"],elements["Γ¹"],kᴹᵛ,fᴹ)
opsh[4](elements["Γ²ₚ"],elements["Γ²"],kᴹᵛ,fᴹ)
opsh[4](elements["Γ³ₚ"],elements["Γ³"],kᴹᵛ,fᴹ)
opsh[4](elements["Γ⁴ₚ"],elements["Γ⁴"],kᴹᵛ,fᴹ)

# d = [zeros(3*nₚ,3*nₚ) zeros(3*nₚ,6*nᵥ) zeros(3*nₚ,9*nᵥ) kᴺᵛ' kᴹᵛ';
#      zeros(6*nᵥ,3*nₚ) kᵋᵋ zeros(6*nᵥ,9*nᵥ) kᴺᵋ' zeros(6*nᵥ,9*nᵥ);
#      zeros(9*nᵥ,3*nₚ) zeros(9*nᵥ,6*nᵥ) kᵏᵏ zeros(9*nᵥ,6*nᵥ) kᴹᵏ';
#      kᴺᵛ kᴺᵋ zeros(6*nᵥ,9*nᵥ) zeros(6*nᵥ,6*nᵥ) zeros(6*nᵥ,9*nᵥ);
#      kᴹᵛ zeros(9*nᵥ,6*nᵥ) kᴹᵏ zeros(9*nᵥ,6*nᵥ) zeros(9*nᵥ,9*nᵥ)]\[-f;zeros(6*nᵥ);zeros(9*nᵥ);fᴺ;fᴹ]

kᵋᵛ = kᴺᵋ\kᴺᵛ
kᵏᵛ = kᴹᵏ\kᴹᵛ
d = (kᵋᵛ'*kᵋᵋ*kᵋᵛ + kᵏᵛ'*kᵏᵏ*kᵏᵛ)\(-f + kᵋᵛ'*kᵋᵋ*(kᴺᵋ\fᴺ) + kᵏᵛ'*kᵏᵏ*(kᴹᵏ\fᴹ))

d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
@save compress=true "jld/patchtest_mix_hr.jld" d₁ d₂ d₃

set∇²𝝭!(elements["Ωᵍ"])
opE = Operator{:Hₑ_ThinShell}(:E=>E,:ν=>ν,:h=>h)
Hₑ, L₂ = opE(elements["Ωᵍ"])

# dᵛ = zeros(3*nₚ)
# for (i,node) in enumerate(nodes)
#     x = Vec{3}((node.x,node.y,node.z))
#     u_ = u(x)
#     dᵛ[3*i-2] = u_[1]
#     dᵛ[3*i-1] = u_[2]
#     dᵛ[3*i]   = u_[3]
# end

# fex = (kᵋᵛ'*kᵋᵋ*kᵋᵛ + kᵏᵛ'*kᵏᵏ*kᵏᵛ)*dᵛ - (-f + kᵋᵛ'*kᵋᵋ*(kᴺᵋ\fᴺ) + kᵏᵛ'*kᵏᵏ*(kᴹᵏ\fᴹ))
# fex1 = (kᵋᵛ'*kᵋᵋ*kᵋᵛ)*dᵛ - (-f + kᵋᵛ'*kᵋᵋ*(kᴺᵋ\fᴺ))
# fex2 = (kᵏᵛ'*kᵏᵏ*kᵏᵛ)*dᵛ - (kᵏᵛ'*kᵏᵏ*(kᴹᵏ\fᴹ))
# fex2 = (kᵏᵛ'*kᵏᵏ*kᵏᵛ)*dᵛ

# dᵋ = zeros(6*nᵥ)
# for i in 1:length(elements["Ω"])
#     dᵋ[18*i-17] = 1.0
#     dᵋ[18*i-14] = 2.0
# end

# dᴺ = zeros(6*nᵥ)
# for i in 1:length(elements["Ω"])
#     dᴺ[18*i-17] = E*h
#     dᴺ[18*i-16] = E*h
#     dᴺ[18*i-14] = E*h
# end

# dᵏ = zeros(9*nᵥ)
# for i in 1:length(elements["Ω"])
#     dᵏ[27*i-24] = -2.
#     dᵏ[27*i-21] = -2.
#     dᵏ[27*i-18] = -4.
#     dᵏ[27*i-24] = 1.
# end

# dᴹ = zeros(9*nᵥ)
# for i in 1:length(elements["Ω"])
#     dᴹ[27*i-24] = -E*h^3/12*2
#     dᴹ[27*i-21] = -E*h^3/12*2
#     dᴹ[27*i-18] = -E*h^3/12*2
#     dᴹ[27*i-24] = E*h^3/12
# end
# println(norm(kᴹᵛ'*dᴹ))
# println(norm(kᵏᵏ*dᵏ+kᴹᵏ'*dᴹ))
# println(norm(kᴹᵛ*dᵛ+kᴹᵏ*dᵏ-fᴹ))
# f1 = kᴹᵛ*dᵛ
# f2 = kᴹᵏ*dᵏ

# println(norm(kᴺᵛ'*dᴺ+kᴹᵛ'*dᴹ))
# println(norm(kᵋᵋ*dᵋ+kᴺᵋ'*dᴺ))
# println(norm(kᵏᵏ*dᵏ+kᴹᵏ'*dᴹ))
# println(norm(kᴺᵛ*dᵛ+kᴺᵋ*dᵋ-fᴺ))
# println(norm(kᴹᵛ*dᵛ+kᴹᵏ*dᵏ-fᴹ))
