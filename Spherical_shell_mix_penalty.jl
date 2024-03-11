using ApproxOperator, JLD

import BenchmarkExample: BenchmarkExample, XLSX
include("import_prescrible_ops.jl")
include("import_Spherical_shell.jl")
ndiv = 32
elements, nodes = import_spherical_mix("msh/sphericalshell_"*string(ndiv)*".msh");

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


d₁ = zeros(nₚ)
d₂ = zeros(nₚ)
d₃ = zeros(nₚ)
d = zeros(3*nₚ)
kᵅ = zeros(3*nₚ,3*nₚ)
fᵅ = zeros(3*nₚ)

kᴳ = zeros(3*nₚ)
ξ = elements["𝐴"][1].𝓖[1]
𝓒 = elements["𝐴"][1].𝓒
N = ξ[:𝝭]
for (i,xᵢ) in enumerate(𝓒)
    I = xᵢ.𝐼
    kᴳ[3*I] = -N[i]*E
end

kᵋᵛ = kᴺᵋ\kᴺᵛ
kᵏᵛ = kᴹᵏ\kᴹᵛ
k = kᵋᵛ'*kᵋᵋ*kᵋᵛ + kᵏᵛ'*kᵏᵏ*kᵏᵛ
op𝐴 = Operator{:SphericalShell_𝐴}()

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
# index = [8,16,24,32]
# for (i,αᵥ) in enumerate([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16])
#     for (j,αᵣ) in enumerate([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16])
        αᵥ = 1e9
        αᵣ = 1e7
        opsα = [
            Operator{:∫vᵢgᵢdΓ}(:α=>αᵥ),
            Operator{:∫δθθdΓ}(:α=>αᵣ),
        ]
        fill!(kᵅ,0.0)
        fill!(fᵅ,0.0)
        opsα[1](elements["Γʳ"],kᵅ,fᵅ)
        opsα[1](elements["Γˡ"],kᵅ,fᵅ)
        opsα[2](elements["Γʳ"],kᵅ,fᵅ)
        opsα[2](elements["Γˡ"],kᵅ,fᵅ)

        d = [k - kᵅ kᴳ;kᴳ' 0]\[-f - fᵅ;0]

        d₁ .= d[1:3:3*nₚ]
        d₂ .= d[2:3:3*nₚ]
        d₃ .= d[3:3:3*nₚ]

        w = op𝐴(elements["𝐴"])
        # println(17*(i-1)+j)
        println(w)

        # XLSX.openxlsx("./xlsx/spherical_penalty_alpha.xlsx", mode="rw") do xf
        #     ind = findfirst((x)->x==ndiv,index)
        #     𝐿₂_row = Char(64+j)*string(i)
        #     xf[ind][𝐿₂_row] = w
        # end
#     end
# end