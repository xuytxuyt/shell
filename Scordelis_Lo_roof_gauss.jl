using ApproxOperator

import BenchmarkExample: BenchmarkExample

include("import_Scordelis_Lo_roof.jl")
ndiv = 51
elements, nodes = import_roof_gauss("msh/scordelislo_"*string(ndiv)*".msh");
nₚ = length(nodes)
s = 2.5*25/ndiv*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

𝑅 = BenchmarkExample.ScordelisLoRoof.𝑅
b₃ = BenchmarkExample.ScordelisLoRoof.𝑞
E = BenchmarkExample.ScordelisLoRoof.𝐸
ν = BenchmarkExample.ScordelisLoRoof.𝜈
h = BenchmarkExample.ScordelisLoRoof.ℎ
cs = BenchmarkExample.cylindricalCoordinate(𝑅)

set𝝭!(elements["Ω"])
set∇𝝭!(elements["Ω"])
set∇²𝝭!(elements["Ω"])
set𝝭!(elements["Ω"])
set𝝭!(elements["Γᵇ"])
set𝝭!(elements["Γʳ"])
set𝝭!(elements["Γᵗ"])
set𝝭!(elements["Γˡ"])
set𝝭!(elements["𝐴"])

prescribe!(elements["Ω"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(ξ¹,ξ²))
prescribe!(elements["Ω"],:a²²=>(ξ¹,ξ²,ξ³)->cs.a²²(ξ¹,ξ²))
prescribe!(elements["Ω"],:a¹²=>(ξ¹,ξ²,ξ³)->cs.a¹²(ξ¹,ξ²))
prescribe!(elements["Ω"],:a₁₁=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(ξ¹,ξ²)[1])
prescribe!(elements["Ω"],:a₁₂=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(ξ¹,ξ²)[2])
prescribe!(elements["Ω"],:a₁₃=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(ξ¹,ξ²)[3])
prescribe!(elements["Ω"],:a₂₁=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(ξ¹,ξ²)[1])
prescribe!(elements["Ω"],:a₂₂=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(ξ¹,ξ²)[2])
prescribe!(elements["Ω"],:a₂₃=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(ξ¹,ξ²)[3])
prescribe!(elements["Ω"],:a₃₁=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(ξ¹,ξ²)[1])
prescribe!(elements["Ω"],:a₃₂=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(ξ¹,ξ²)[2])
prescribe!(elements["Ω"],:a₃₃=>(ξ¹,ξ²,ξ³)->cs.𝒂₃(ξ¹,ξ²)[3])
prescribe!(elements["Ω"],:Γ¹₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₁(ξ¹,ξ²))
prescribe!(elements["Ω"],:Γ¹₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₁₂(ξ¹,ξ²))
prescribe!(elements["Ω"],:Γ¹₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ¹₂₂(ξ¹,ξ²))
prescribe!(elements["Ω"],:Γ²₁₁=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₁(ξ¹,ξ²))
prescribe!(elements["Ω"],:Γ²₁₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₁₂(ξ¹,ξ²))
prescribe!(elements["Ω"],:Γ²₂₂=>(ξ¹,ξ²,ξ³)->cs.Γ²₂₂(ξ¹,ξ²))
prescribe!(elements["Ω"],:b₁=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Ω"],:b₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Ω"],:b₃=>(ξ¹,ξ²,ξ³)->b₃)
prescribe!(elements["Γᵇ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:n₁₁=>(ξ¹,ξ²,ξ³)->1.0)
prescribe!(elements["Γᵇ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:n₂₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵇ"],:n₃₃=>(ξ¹,ξ²,ξ³)->1.0)
prescribe!(elements["Γᵗ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:n₁₁=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:n₂₂=>(ξ¹,ξ²,ξ³)->1.0)
prescribe!(elements["Γᵗ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γᵗ"],:n₃₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:g₁=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:g₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:g₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:n₁₁=>(ξ¹,ξ²,ξ³)->1.0)
prescribe!(elements["Γˡ"],:n₁₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:n₁₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:n₂₂=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:n₂₃=>(ξ¹,ξ²,ξ³)->0.0)
prescribe!(elements["Γˡ"],:n₃₃=>(ξ¹,ξ²,ξ³)->0.0)
ops = [
    Operator{:∫εᵢⱼNᵢⱼκᵢⱼMᵢⱼdΩ}(:E=>E,:ν=>ν,:h=>h,),
    Operator{:∫vᵢbᵢdΩ}(),
    Operator{:∫vᵢgᵢdΓ}(:α=>1e9*E),
    Operator{:ScordelisLoRoof_𝐴}()
]
k = zeros(3*nₚ,3*nₚ)
f = zeros(3*nₚ)

ops[1](elements["Ω"],k)
ops[2](elements["Ω"],f)
ops[3](elements["Γᵇ"],k,f)
ops[3](elements["Γᵗ"],k,f)
ops[3](elements["Γˡ"],k,f)

d = k\f
d₁ = d[1:3:3*nₚ]
d₂ = d[2:3:3*nₚ]
d₃ = d[3:3:3*nₚ]

push!(nodes,:d₁=>d₁,:d₂=>d₂,:d₃=>d₃)
w = ops[4](elements["𝐴"])
