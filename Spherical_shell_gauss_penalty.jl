using ApproxOperator, JLD

import BenchmarkExample: BenchmarkExample

include("import_Spherical_shell.jl")

𝑅 = BenchmarkExample.SphericalShell.𝑅
F = BenchmarkExample.SphericalShell.𝐹
E = BenchmarkExample.SphericalShell.𝐸
ν = BenchmarkExample.SphericalShell.𝜈
h = BenchmarkExample.SphericalShell.ℎ
𝜃 =  BenchmarkExample.SphericalShell.𝜃₂
cs = BenchmarkExample.sphericalCoordinate(𝑅)

ndiv = 64
elements, nodes = import_Spherical_shell_gauss("msh/sphericalshell_"*string(ndiv)*".msh");
nₚ = length(nodes)
s = 2.5*𝑅*𝜃/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set∇²𝝭!(elements["Ω"])
set𝝭!(elements["Γᵇ"])
set∇𝝭!(elements["Γʳ"])
set𝝭!(elements["Γᵗ"])
set∇𝝭!(elements["Γˡ"])
set𝝭!(elements["𝐴"])
set𝝭!(elements["𝐵"])

eval(prescribleBoundary)

ops = [
    Operator{:∫εᵢⱼNᵢⱼκᵢⱼMᵢⱼdΩ}(:E=>E,:ν=>ν,:h=>h),
    Operator{:∫vᵢtᵢdΓ}(),
    Operator{:SphericalShell_𝐴}(),
    Operator{:SphericalShell_𝐵}(),
]

k = zeros(3*nₚ,3*nₚ)
f = zeros(3*nₚ)

ops[1](elements["Ω"],k)
ops[2](elements["𝐴"],f)
ops[2](elements["𝐵"],f)

# for (i,αᵥ) in enumerate([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16])
#     for (j,αᵣ) in enumerate([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16])
# for (i,αᵥ) in enumerate([1e9])
#     for (j,αᵣ) in enumerate([1e5])
αᵥ = 1e9
αᵣ = 1e7
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

        d = (k+kᵅ)\(f+fᵅ)
        d₁ = d[1:3:3*nₚ]
        d₂ = d[2:3:3*nₚ]
        d₃ = d[3:3:3*nₚ]

        push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
        w₁ = ops[3](elements["𝐴"])
        w₂ = ops[4](elements["𝐵"])

        println(w₁)
        println(w₂)
        # @save compress=true "jld/Spherical_shell_gauss_penalty_"*string(ndiv)*".jld" d₁ d₂ d₃
#     end
# end
