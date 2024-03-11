using ApproxOperator, JLD, TimerOutputs
const to = TimerOutput()
import BenchmarkExample: BenchmarkExample
include("import_prescrible_ops.jl")
include("import_Spherical_shell.jl")

ndiv = 8
nₘ = 55
𝑅 = BenchmarkExample.SphericalShell.𝑅
E = BenchmarkExample.SphericalShell.𝐸
ν = BenchmarkExample.SphericalShell.𝜈
h = BenchmarkExample.SphericalShell.ℎ
𝜃 =  BenchmarkExample.SphericalShell.𝜃₂
𝐹 = BenchmarkExample.SphericalShell.𝐹

cs = BenchmarkExample.sphericalCoordinate(𝑅)

elements, nodes = import_spherical_gauss("msh/sphericalshell_"*string(ndiv)*".msh");
nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])*3)
# s = 2.26*𝑅*𝜃/(ndiv-1)*ones(nₚ)
s = 2.5*𝑅*𝜃/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)
eval(prescribeForGauss)
@timeit to "shape gauss" begin
set∇²𝝭!(elements["Ω"])
end
eval(opsGauss)
k = zeros(3*nₚ,3*nₚ)
@timeit to "Gauss" begin
op(elements["Ω"],k)
end
elements, nodes = import_spherical_mix("msh/sphericalshell_"*string(ndiv)*".msh");

nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])*3)
# s = 2.26*𝑅*𝜃/(ndiv-1)*ones(nₚ)
s = 2.5*𝑅*𝜃/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

eval(prescribeForMix)
eval(prescribeForNitsche)
eval(prescribeForPenalty)
eval(prescribeVariables)

set𝝭!(elements["Ω"])
set∇²𝝭!(elements["Ωₚ"])
set∇𝝭!(elements["Γₚ"])
set∇𝝭!(elements["Ωₘ"])
@timeit to "shape mix" begin
set∇𝝭!(elements["Γₘ"])
end

@timeit to "shape HR" begin
set∇𝝭!(elements["Γʳ"])
set∇𝝭!(elements["Γˡ"])
end

@timeit to "shape nitsche" begin
set∇̂³𝝭!(elements["Γʳ"])
set∇̂³𝝭!(elements["Γˡ"])
set∇𝝭!(elements["Γʳ"])
set∇𝝭!(elements["Γˡ"])
end

@timeit to "shape penalty" begin
set∇𝝭!(elements["Γʳ"])
set∇𝝭!(elements["Γˡ"])
end


kᵋᵋ = zeros(6*nᵥ,6*nᵥ)
kᴺᵋ = zeros(6*nᵥ,6*nᵥ)
kᴺᵛ = zeros(6*nᵥ,3*nₚ)
kᵏᵏ = zeros(9*nᵥ,9*nᵥ)
kᴹᵏ = zeros(9*nᵥ,9*nᵥ)
kᴹᵛ = zeros(9*nᵥ,3*nₚ)

eval(opsMix)
@timeit to "RKGSI" begin
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

end
kᵋᵛ = kᴺᵋ\kᴺᵛ
kᵏᵛ = kᴹᵏ\kᴹᵛ
k = kᵋᵛ'*kᵋᵋ*kᵋᵛ + kᵏᵛ'*kᵏᵏ*kᵏᵛ

eval(opsHR)
fᴺ = zeros(6*nᵥ)
fᴹ = zeros(9*nᵥ)
@timeit to "HR" begin
opsh[1](elements["Γʳₚ"],elements["Γʳ"],kᴺᵛ,fᴺ)
opsh[1](elements["Γˡₚ"],elements["Γˡ"],kᴺᵛ,fᴺ)
opsh[2](elements["Γʳₚ"],elements["Γʳ"],kᴹᵛ,fᴹ)
opsh[2](elements["Γˡₚ"],elements["Γˡ"],kᴹᵛ,fᴹ)
opsh[3](elements["Γʳₚ"],elements["Γʳ"],kᴹᵛ,fᴹ)
opsh[3](elements["Γˡₚ"],elements["Γˡ"],kᴹᵛ,fᴹ)
kᵋᵛ'*kᵋᵋ*(kᴺᵋ\fᴺ) + kᵏᵛ'*kᵏᵏ*(kᴹᵏ\fᴹ)
end

eval(opsNitsche)
kᵛ = zeros(3*nₚ,3*nₚ)
fᵛ = zeros(3*nₚ)
@timeit to "Nitsche" begin
opsv[1](elements["Γʳ"],kᵛ,fᵛ)
opsv[1](elements["Γˡ"],kᵛ,fᵛ)
opsv[2](elements["Γʳ"],kᵛ,fᵛ)
opsv[2](elements["Γˡ"],kᵛ,fᵛ)
opsv[3](elements["Γʳ"],kᵛ,fᵛ)
opsv[3](elements["Γˡ"],kᵛ,fᵛ)
end

αᵥ = 1.0
αᵣ = 1.0
eval(opsPenalty)
kᵅ = zeros(3*nₚ,3*nₚ)
fᵅ = zeros(3*nₚ)
@timeit to "Penalty" begin
opsα[1](elements["Γʳ"],kᵅ,fᵅ)
opsα[1](elements["Γˡ"],kᵅ,fᵅ)
opsα[2](elements["Γʳ"],kᵅ,fᵅ)
opsα[2](elements["Γˡ"],kᵅ,fᵅ)
end
show(to)