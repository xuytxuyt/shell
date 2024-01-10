using ApproxOperator

import BenchmarkExample: BenchmarkExample

include("import_Scordelis_Lo_roof.jl")
ndiv = 10
elements, nodes = import_roof_gauss("msh/scordelislo_"*string(ndiv)*".msh");

𝑅 = BenchmarkExample.ScordelisLoRoof.𝑅
b₃ = BenchmarkExample.ScordelisLoRoof.𝑞
cs = BenchmarkExample.cylindricalCoordinate(𝑅)
prescribe!(elements["Ω"],:a¹¹=>(ξ¹,ξ²,ξ³)->cs.a¹¹(ξ¹,ξ²))
prescribe!(elements["Ω"],:a₁₁=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(ξ¹,ξ²)[1])
prescribe!(elements["Ω"],:a₁₂=>(ξ¹,ξ²,ξ³)->cs.𝒂₁(ξ¹,ξ²)[1])
prescribe!(elements["Ω"],:a₂₁=>(ξ¹,ξ²,ξ³)->cs.𝒂₂(ξ¹,ξ²)[1])

ops = [
    Operator{:∫εᵢⱼNᵢⱼκᵢⱼMᵢⱼdΩ}
    Operator{:∫vᵢbᵢdΩ}
    Operator{:∫vᵢgᵢdΓ}
]