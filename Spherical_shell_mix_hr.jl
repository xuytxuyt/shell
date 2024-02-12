using ApproxOperator, JLD

import BenchmarkExample: BenchmarkExample
include("import_prescrible_ops.jl")
include("import_Spherical_shell.jl")
ndiv = 8
elements, nodes = import_roof_mix("msh/sphericalshell_"*string(ndiv)*".msh");

nₘ = 21
𝑅 = BenchmarkExample.SphericalShell.𝑅
E = BenchmarkExample.SphericalShell.𝐸
ν = BenchmarkExample.SphericalShell.𝜈
h = BenchmarkExample.SphericalShell.ℎ
𝜃 =  BenchmarkExample.SphericalShell.𝜃₂
𝐹 = BenchmarkExample.SphericalShell.𝐹

cs = BenchmarkExample.sphericalCoordinate(𝑅)
nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])*3)
s = 2.5*𝑅*𝜃/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

eval(prescribeForMix)
eval(prescribeVariables)

set𝝭!(elements["Ω"])
set∇²𝝭!(elements["Ωₚ"])
set∇𝝭!(elements["Γₚ"])
set∇𝝭!(elements["Ωₘ"])
set∇𝝭!(elements["Γₘ"])
set∇𝝭!(elements["Γʳ"])
set∇𝝭!(elements["Γˡ"])
set𝝭!(elements["𝐴"])
set𝝭!(elements["𝐵"])

opForce = Operator{:∫vᵢtᵢdΓ}()
f = zeros(3*nₚ)
opForce(elements["𝐴"],f)
opForce(elements["𝐵"],f)

eval(opsMix)

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
opsh[1](elements["Γʳₚ"],elements["Γʳ"],kᴺᵛ,fᴺ)
opsh[1](elements["Γˡₚ"],elements["Γˡ"],kᴺᵛ,fᴺ)
opsh[2](elements["Γʳₚ"],elements["Γʳ"],kᴹᵛ,fᴹ)
opsh[2](elements["Γˡₚ"],elements["Γˡ"],kᴹᵛ,fᴹ)
opsh[3](elements["Γʳₚ"],elements["Γʳ"],kᴹᵛ,fᴹ)
opsh[3](elements["Γˡₚ"],elements["Γˡ"],kᴹᵛ,fᴹ)

opα = Operator{:∫vᵢgᵢdΓ}(:α=>1e4*E)
kᵅ = zeros(3*nₚ,3*nₚ)
fᵅ = zeros(3*nₚ)
opα(elements["𝐴"],kᵅ,fᵅ)

kᵋᵛ = kᴺᵋ\kᴺᵛ
kᵏᵛ = kᴹᵏ\kᴹᵛ
d = (kᵋᵛ'*kᵋᵋ*kᵋᵛ + kᵏᵛ'*kᵏᵏ*kᵏᵛ)\(-f + kᵋᵛ'*kᵋᵋ*(kᴺᵋ\fᴺ) + kᵏᵛ'*kᵏᵏ*(kᴹᵏ\fᴹ))

d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

op𝐴 = Operator{:SphericalShell_𝐴}()
push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
w = op𝐴(elements["𝐴"])

println(w)
@save compress=true "jld/spherical_shell_mix_hr_"*string(ndiv)*".jld" d₁ d₂ d₃