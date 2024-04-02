using ApproxOperator, JLD

import BenchmarkExample: BenchmarkExample
include("import_prescrible_ops.jl")
include("import_Spherical_shell.jl")

nₘ = 21
𝑅 = BenchmarkExample.SphericalShell.𝑅
F = BenchmarkExample.SphericalShell.𝐹
E = BenchmarkExample.SphericalShell.𝐸
ν = BenchmarkExample.SphericalShell.𝜈
h = BenchmarkExample.SphericalShell.ℎ
𝜃 =  BenchmarkExample.SphericalShell.𝜃₂
cs = BenchmarkExample.sphericalCoordinate(𝑅)

ndiv = 32
elements, nodes = import_spherical_gauss("msh/sphericalshell_"*string(ndiv)*".msh");
nₚ = length(nodes)
s = 2.5*𝑅*𝜃/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

eval(prescribeForGauss)
eval(prescribeForNitsche)
eval(prescribeVariables)
eval(opsGauss)

set∇²𝝭!(elements["Ω"])
set∇𝝭!(elements["Γʳ"])
set∇𝝭!(elements["Γˡ"])
set𝝭!(elements["𝐴"])
set𝝭!(elements["𝐵"])


opForce = Operator{:∫vᵢtᵢdΓ}()
f = zeros(3*nₚ)
opForce(elements["𝐴"],f)
opForce(elements["𝐵"],f)

k = zeros(3*nₚ,3*nₚ)

op(elements["Ω"],k)

kᴳ = zeros(3*nₚ)
ξ = elements["𝐴"][1].𝓖[1]
𝓒 = elements["𝐴"][1].𝓒
N = ξ[:𝝭]
for (i,xᵢ) in enumerate(𝓒)
    I = xᵢ.𝐼
    kᴳ[3*I] = -N[i]*1e5*E
end

eval(opsNitsche)
# kᵛ = zeros(3*nₚ,3*nₚ)
# fᵛ = zeros(3*nₚ)
opsv[1](elements["Γʳ"],k,f)
opsv[1](elements["Γˡ"],k,f)
opsv[2](elements["Γʳ"],k,f)
opsv[2](elements["Γˡ"],k,f)
opsv[3](elements["Γʳ"],k,f)
opsv[3](elements["Γˡ"],k,f)

# for (i,αᵥ) in enumerate([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16])
#     for (j,αᵣ) in enumerate([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16])
# for (i,αᵥ) in enumerate([1e9])
#     for (j,αᵣ) in enumerate([1e5])
αᵥ = 1e5
αᵣ = 1e3
        opΓ = [
            Operator{:∫vᵢgᵢdΓ}(:α=>αᵥ*E),
            Operator{:∫δθθdΓ}(:α=>αᵣ*E)
        ]
        kᵅ = zeros(3*nₚ,3*nₚ)
        fᵅ = zeros(3*nₚ)
        opΓ[1](elements["Γʳ"],kᵅ,fᵅ)
        opΓ[1](elements["Γˡ"],kᵅ,fᵅ)
        opΓ[2](elements["Γʳ"],kᵅ,fᵅ)
        opΓ[2](elements["Γˡ"],kᵅ,fᵅ)

        d = [k+kᵅ kᴳ;kᴳ' 0]\[f+fᵅ;0]
        d₁ = d[1:3:3*nₚ]
        d₂ = d[2:3:3*nₚ]
        d₃ = d[3:3:3*nₚ]

        op𝐴 = Operator{:SphericalShell_𝐴}()
        push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
        w = op𝐴(elements["𝐴"])

        # @save compress=true "jld/Spherical_shell_gauss_penalty_"*string(ndiv)*".jld" d₁ d₂ d₃
#     end
# end
