using ApproxOperator, JLD, XLSX 

import BenchmarkExample: BenchmarkExample
include("import_prescrible_ops.jl")
include("import_Scordelis_Lo_roof.jl")
ndiv = 24
elements, nodes = import_roof_mix("msh/scordelislo_"*string(ndiv)*".msh");

nₘ = 55
𝑅 = BenchmarkExample.ScordelisLoRoof.𝑅
𝐿 = BenchmarkExample.ScordelisLoRoof.𝐿
b₃ = BenchmarkExample.ScordelisLoRoof.𝑞
E = BenchmarkExample.ScordelisLoRoof.𝐸
ν = BenchmarkExample.ScordelisLoRoof.𝜈
h = BenchmarkExample.ScordelisLoRoof.ℎ
cs = BenchmarkExample.cylindricalCoordinate(𝑅)
nₚ = length(nodes)
nᵥ = Int(length(elements["Ω"])*3)
s = 2.5*𝐿/2/(ndiv-1)*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

eval(prescribeForMix)
eval(prescribeVariables)
eval(prescribeForNitsche)

set𝝭!(elements["Ω"])
set∇²𝝭!(elements["Ωₚ"])
set∇𝝭!(elements["Γₚ"])
set∇𝝭!(elements["Ωₘ"])
set∇𝝭!(elements["Γₘ"])
set∇̂³𝝭!(elements["Γᵇ"])
set∇̂³𝝭!(elements["Γᵗ"])
set∇̂³𝝭!(elements["Γˡ"])
set∇𝝭!(elements["Γᵇ"])
set∇𝝭!(elements["Γᵗ"])
set∇𝝭!(elements["Γˡ"])
set𝝭!(elements["𝐴"])

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
fᴺ = zeros(6*nᵥ)
fᴹ = zeros(3*nᵥ)

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

kᵋᵛ = kᴺᵋ\kᴺᵛ
kᵏᵛ = kᴹᵏ\kᴹᵛ
op𝐴 = Operator{:ScordelisLoRoof_𝐴}()

eval(opsNitsche)
kᵛ = zeros(3*nₚ,3*nₚ)
fᵛ = zeros(3*nₚ)
opsv[1](elements["Γᵇ"],kᵛ,fᵛ)
opsv[1](elements["Γᵗ"],kᵛ,fᵛ)
opsv[1](elements["Γˡ"],kᵛ,fᵛ)
opsv[2](elements["Γᵇ"],kᵛ,fᵛ)
opsv[2](elements["Γᵗ"],kᵛ,fᵛ)
opsv[2](elements["Γˡ"],kᵛ,fᵛ)
opsv[3](elements["Γᵗ"],kᵛ,fᵛ)
opsv[3](elements["Γˡ"],kᵛ,fᵛ)

d₁ = zeros(nₚ)
d₂ = zeros(nₚ)
d₃ = zeros(nₚ)
d = zeros(3*nₚ)
kᵅ = zeros(3*nₚ,3*nₚ)
fᵅ = zeros(3*nₚ)
k_ = kᵋᵛ'*kᵋᵋ*kᵋᵛ + kᵏᵛ'*kᵏᵏ*kᵏᵛ - kᵛ
f_ = -f - fᵛ
push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
index = [8,16,24,32]
for (i,αᵥ) in enumerate([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16])
    for (j,αᵣ) in enumerate([1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9,1e10,1e11,1e12,1e13,1e14,1e15,1e16])
        opsα = [
            Operator{:∫vᵢgᵢdΓ}(:α=>αᵥ*E),
            Operator{:∫δθθdΓ}(:α=>αᵣ*E),
        ]
        fill!(kᵅ,0.0)
        fill!(fᵅ,0.0)
        opsα[1](elements["Γᵇ"],kᵅ,fᵅ)
        opsα[1](elements["Γᵗ"],kᵅ,fᵅ)
        opsα[1](elements["Γˡ"],kᵅ,fᵅ)
        opsα[2](elements["Γᵗ"],kᵅ,fᵅ)
        opsα[2](elements["Γˡ"],kᵅ,fᵅ)

        d .= (k_ - kᵅ)\(f_ - fᵅ)

        d₁ .= d[1:3:3*nₚ]
        d₂ .= d[2:3:3*nₚ]
        d₃ .= d[3:3:3*nₚ]

        w = op𝐴(elements["𝐴"])
        println(17*(i-1)+j)
        println(w)

        XLSX.openxlsx("./xlsx/scordelis_lo_roof_nitsche_alpha.xlsx", mode="rw") do xf
            ind = findfirst((x)->x==ndiv,index)
            𝐿₂_row = Char(64+j)*string(i)
            xf[ind][𝐿₂_row] = w
        end
    end
end