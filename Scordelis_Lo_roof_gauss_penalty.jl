using ApproxOperator, JLD

import BenchmarkExample: BenchmarkExample

include("import_Scordelis_Lo_roof.jl")

𝑅 = BenchmarkExample.ScordelisLoRoof.𝑅
𝐿 = BenchmarkExample.ScordelisLoRoof.𝐿
b₃ = BenchmarkExample.ScordelisLoRoof.𝑞
E = BenchmarkExample.ScordelisLoRoof.𝐸
ν = BenchmarkExample.ScordelisLoRoof.𝜈
h = BenchmarkExample.ScordelisLoRoof.ℎ
cs = BenchmarkExample.cylindricalCoordinate(𝑅)

ndiv = 16
elements, nodes = import_roof_gauss("msh/scordelislo_"*string(ndiv)*".msh");
nₚ = length(nodes)
s = 3.5*𝐿/2/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set∇²𝝭!(elements["Ω"])
set𝝭!(elements["Γᵇ"])
set𝝭!(elements["Γʳ"])
set∇𝝭!(elements["Γᵗ"])
set∇𝝭!(elements["Γˡ"])
set𝝭!(elements["𝐴"])

eval(prescribleBoundary)

ops = [
    Operator{:∫εᵢⱼNᵢⱼκᵢⱼMᵢⱼdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫vᵢbᵢdΩ}(),
    Operator{:ScordelisLoRoof_𝐴}()
]

k = zeros(3*nₚ,3*nₚ)
f = zeros(3*nₚ)

ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)

for (i,αᵥ) in enumerate([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16])
    for (j,αᵣ) in enumerate([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16])
        opΓ = [
            Operator{:∫vᵢgᵢdΓ}(:α=>αᵥ*E),
            Operator{:∫δθθdΓ}(:α=>αᵣ*E)
        ]
        kᵅ = zeros(3*nₚ,3*nₚ)
        fᵅ = zeros(3*nₚ)
        ops[3](elements["Γᵇ"],kᵅ,fᵅ)
        ops[3](elements["Γᵗ"],kᵅ,fᵅ)
        ops[3](elements["Γˡ"],kᵅ,fᵅ)
        ops[4](elements["Γᵗ"],kᵅ,fᵅ)
        ops[4](elements["Γˡ"],kᵅ,fᵅ)

        d = (k+kᵅ)\(f+fᵅ)
        d₁ = d[1:3:3*nₚ]
        d₂ = d[2:3:3*nₚ]
        d₃ = d[3:3:3*nₚ]

        push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
        w = ops[5](elements["𝐴"])

        println(w)
    end
end

# @save compress=true "jld/scordelislo_gauss_"*string(ndiv)*".jld" d₁ d₂ d₃