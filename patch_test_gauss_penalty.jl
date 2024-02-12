
using ApproxOperator, Tensors, JLD

import BenchmarkExample: BenchmarkExample
include("import_prescrible_ops.jl")
include("import_patch_test.jl")
elements, nodes = import_roof_gauss("msh/patchtest.msh")

𝑅 = BenchmarkExample.PatchTestThinShell.𝑅
𝐿 = BenchmarkExample.PatchTestThinShell.𝐿
h = BenchmarkExample.PatchTestThinShell.ℎ
E = BenchmarkExample.PatchTestThinShell.𝐸
ν = BenchmarkExample.PatchTestThinShell.𝜈

nₘ = 55
cs = BenchmarkExample.cylindricalCoordinate(𝑅)
# cs = BenchmarkExample.cartesianCoordinate()

nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])*3)
s = 2.5*0.1*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

n = 2
u(x) = Vec{3}(((x[1]+2*x[2])^n,(3*x[1]+4*x[2])^n,(5*x[1]+6*x[2])^n))
vs = BenchmarkExample.PatchTestThinShell.variables(cs,u)

eval(prescribeForGauss)
eval(prescribeForPenalty)
eval(prescribeVariables)
eval(opsGauss)

set∇²𝝭!(elements["Ω"])
set∇𝝭!(elements["Γ¹"])
set∇𝝭!(elements["Γ²"])
set∇𝝭!(elements["Γ³"])
set∇𝝭!(elements["Γ⁴"])

opForce = Operator{:∫vᵢbᵢdΩ}()
f = zeros(3*nₚ)
opForce(elements["Ω"],f)

k = zeros(3*nₚ,3*nₚ)
op(elements["Ω"],k)

αᵥ = 1e9
αᵣ = 1e9
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

d = (k + kᵅ)\(f + fᵅ)

d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
@save compress=true "jld/patchtest_gauss_penalty.jld" d₁ d₂ d₃

set∇²𝝭!(elements["Ωᵍ"])
opE = Operator{:Hₑ_ThinShell}(:E=>E,:ν=>ν,:h=>h)
Hₑ, L₂ = opE(elements["Ωᵍ"])