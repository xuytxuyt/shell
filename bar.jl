using Revise, ApproxOperator, LinearAlgebra, XLSX

include("import_Scordelis_Lo_roof.jl")

ndiv= 3
elements,nodes = import_roof_bar("./msh/bar_"*string(ndiv)*".msh")
nₚ = length(nodes)

s = 1.5*π/2/ndiv*ones(nₚ)
push!(nodes,:s₁=>s,:s₂=>s,:s₃=>s)

set𝝭!(elements["Ω"])
set𝝭!(elements["Γᵍ"])

i=1/10
R = 1
h = R*i
E = 3e6
I = h^3/12
A = h
EI = E*I
EA = E*A
P  = 1000
a₁ = (sin(s/R),cos(s/R))
a₃ = (-cos(s/R),sin(s/R))

ApproxOperator.prescribe!(elements["Γᵍ"],:g₁=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:g₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:g₃=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₁=>(x,y,z)->1.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:n₁₂=>(x,y,z)->0.0)
ApproxOperator.prescribe!(elements["Γᵍ"],:n₂₂=>(x,y,z)->1.0)

ops = [
       Operator{:∫∫κεds}(:EI=>EI,:EA=>EA,:a₁=>a₁,:a₃=>a₃),
       Operator{:∫∫vᵢgᵢdsds}(:α=>1e9*E),
]
k = zeros(2*nₚ,2*nₚ)
f = zeros(2*nₚ)

ops[1](elements["Ω"],k)
ops[2](elements["Γᵍ"],k,f)
d = k\f

