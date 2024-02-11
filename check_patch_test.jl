using Tensors
import BenchmarkExample: BenchmarkExample

# cs = BenchmarkExample.cylindricalCoordinate(1.0)
cs = BenchmarkExample.cartesianCoordinate()

BenchmarkExample.PatchTestThinShell.𝜈 = 0.0

n = 2
# u(x) = Vec{3}(((x[1]+2*x[2])^n,(3*x[1]+4*x[2])^n,(5*x[1]+6*x[2])^n))
# u(x) = Vec{3}((x[1]^n,0.0,(5*x[1]+6*x[2])^n))
# u(x) = Vec{3}((x[1]^n,0.0,0.0))
u(x) = Vec{3}((0.0,0.0,(x[1]+x[2])^n))

vs = BenchmarkExample.PatchTestThinShell.variables(cs,u)

x = Vec{2}((rand(),rand()))

println(vs.𝜃ₙ(x,rand(),rand()))
println(vs.𝜺(x))
println(vs.𝜿(x))
println(vs.𝑪(x))
println(vs.𝑵(x))
println(vs.𝑴(x))
println(vs.𝒃(x))
∇𝑴 = gradient(vs.𝑴,x)
div𝑴 = ∇𝑴[:,1,1]+∇𝑴[:,2,2]
println(∇𝑴)
println(div𝑴)
# println(vs.∇𝒖(x)[:,1])
# println(vs.∇𝒖(x)[:,2])