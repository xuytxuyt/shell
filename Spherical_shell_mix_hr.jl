using ApproxOperator, JLD, TimerOutputs
const to = TimerOutput()
import BenchmarkExample: BenchmarkExample
include("import_prescrible_ops.jl")
include("import_Spherical_shell.jl")
ndiv = 32
elements, nodes = import_spherical_mix("msh/sphericalshell_"*string(ndiv)*".msh");

nₘ = 21
𝑅 = BenchmarkExample.SphericalShell.𝑅
E = BenchmarkExample.SphericalShell.𝐸
ν = BenchmarkExample.SphericalShell.𝜈
h = BenchmarkExample.SphericalShell.ℎ
𝜃₁ =  BenchmarkExample.SphericalShell.𝜃₂
𝜃₂ =  BenchmarkExample.SphericalShell.𝜃₁
𝐹 = BenchmarkExample.SphericalShell.𝐹

cs = BenchmarkExample.sphericalCoordinate(𝑅)
nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])*3)
# s = 2.5*𝑅*𝜃₁/(ndiv-1)*ones(nₚ)
# push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
s₁ = 4*𝑅*𝜃₁/(ndiv-1)*ones(nₚ)
s₂ = 4*𝑅*𝜃₂/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s₁,:s₂=>s₂,:s₃=>s₁)

eval(prescribeForMix)
eval(prescribeForPenalty)
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
# f[elements["𝐴"][1].𝓒[1].𝐼] += 1.0
# f[elements["𝐵"][1].𝓒[1].𝐼] -= 1.0

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

# αᵥ = 1e7*E
# αᵣ = 1e3*E
# eval(opsPenalty)
# kᵅ = zeros(3*nₚ,3*nₚ)
# fᵅ = zeros(3*nₚ)
# opsα[1](elements["Γˡ"],kᵅ,fᵅ)
# opsα[1](elements["Γʳ"],kᵅ,fᵅ)
# opsα[2](elements["Γˡ"],kᵅ,fᵅ)
# opsα[2](elements["Γʳ"],kᵅ,fᵅ)

# opsα[1](elements["𝐴"],kᵅ,fᵅ)
kᴳ = zeros(3*nₚ)
ξ = elements["𝐴"][1].𝓖[1]
𝓒 = elements["𝐴"][1].𝓒
N = ξ[:𝝭]
for (i,xᵢ) in enumerate(𝓒)
    I = xᵢ.𝐼
    kᴳ[3*I] = -N[i]*1e0*E
end

kᵋᵛ = kᴺᵋ\kᴺᵛ
kᵏᵛ = kᴹᵏ\kᴹᵛ
k = kᵋᵛ'*kᵋᵋ*kᵋᵛ + kᵏᵛ'*kᵏᵏ*kᵏᵛ
d = [k kᴳ;kᴳ' 0]\[-f;0]
# d = -k\f

d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

op𝐴 = Operator{:SphericalShell_𝐴}()
push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
w = op𝐴(elements["𝐴"])

println(w)
@save compress=true "jld/spherical_shell_mix_hr_"*string(ndiv)*".jld" d₁ d₂ d₃