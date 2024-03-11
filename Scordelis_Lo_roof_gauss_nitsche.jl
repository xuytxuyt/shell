using ApproxOperator, JLD

import BenchmarkExample: BenchmarkExample

include("import_Scordelis_Lo_roof.jl")
include("import_prescrible_ops.jl")

𝑅 = BenchmarkExample.ScordelisLoRoof.𝑅
𝐿 = BenchmarkExample.ScordelisLoRoof.𝐿
b₃ = BenchmarkExample.ScordelisLoRoof.𝑞
E = BenchmarkExample.ScordelisLoRoof.𝐸
ν = BenchmarkExample.ScordelisLoRoof.𝜈
h = BenchmarkExample.ScordelisLoRoof.ℎ
cs = BenchmarkExample.cylindricalCoordinate(𝑅)

ndiv = 20
elements, nodes = import_roof_gauss("msh/scordelislo_"*string(ndiv)*".msh");
nₚ = length(nodes)
s = 2.5*𝐿/2/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)


eval(prescribeForGauss)
eval(prescribeVariables)
eval(prescribeForNitsche)

set∇²𝝭!(elements["Ω"])
set∇̂³𝝭!(elements["Γᵇ"])
set∇̂³𝝭!(elements["Γᵗ"])
set∇̂³𝝭!(elements["Γˡ"])
set∇𝝭!(elements["Γᵇ"])
set∇𝝭!(elements["Γᵗ"])
set∇𝝭!(elements["Γˡ"])
set𝝭!(elements["𝐴"])
eval(opsGauss)
opForce = Operator{:∫vᵢbᵢdΩ}()
ops𝐴 = Operator{:ScordelisLoRoof_𝐴}()

k = zeros(3*nₚ,3*nₚ)
f = zeros(3*nₚ)

op(elements["Ω"],k)
opForce(elements["Ω"],f)

eval(opsNitsche)
# opsv[1](elements["Γᵇ"],k,f)
# opsv[1](elements["Γᵗ"],k,f)
# opsv[1](elements["Γˡ"],k,f)
# opsv[2](elements["Γᵇ"],k,f)
# opsv[2](elements["Γᵗ"],k,f)
# opsv[2](elements["Γˡ"],k,f)
# opsv[3](elements["Γᵗ"],k,f)
# opsv[3](elements["Γˡ"],k,f)

# for (i,αᵥ) in enumerate([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16])
#     for (j,αᵣ) in enumerate([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16])
# for (i,αᵥ) in enumerate([1e6])
#     for (j,αᵣ) in enumerate([1e6])
αᵥ = 1e7
αᵣ = 1e6
        opΓ = [
            Operator{:∫vᵢgᵢdΓ}(:α=>αᵥ),
            Operator{:∫δθθdΓ}(:α=>αᵣ)
        ]
        kᵅ = zeros(3*nₚ,3*nₚ)
        fᵅ = zeros(3*nₚ)

        opΓ[1](elements["Γᵇ"],kᵅ,fᵅ)
        opΓ[1](elements["Γᵗ"],kᵅ,fᵅ)
        opΓ[1](elements["Γˡ"],kᵅ,fᵅ)
        opΓ[2](elements["Γᵗ"],kᵅ,fᵅ)
        opΓ[2](elements["Γˡ"],kᵅ,fᵅ)

        d = (k+kᵅ)\(f+fᵅ)
        d₁ = d[1:3:3*nₚ]
        d₂ = d[2:3:3*nₚ]
        d₃ = d[3:3:3*nₚ]

        push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
        w = op𝐴(elements["𝐴"])

        println(w)

        @save compress=true "jld/scordelislo_gauss_nitsche_"*string(ndiv)*".jld" d₁ d₂ d₃
#     end
# end