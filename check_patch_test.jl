using Tensors
import BenchmarkExample: BenchmarkExample

cs = BenchmarkExample.cylindricalCoordinate(1.0)

BenchmarkExample.PatchTestThinShell.𝜈 = 0.0

n = 1
u(x) = Vec{3}(((x[1]+2*x[2])^n,(3*x[1]+4*x[2])^n,(5*x[1]+6*x[2])^n))

vs = BenchmarkExample.PatchTestThinShell.variables(cs,u)

x = Vec{2}((rand(),rand()))

println(vs.𝜃ₙ(x,rand(),rand()))
println(vs.𝜺(x))
println(vs.𝜿(x))
println(vs.𝑵(x))
println(vs.𝑴(x))
println(vs.𝒃(x))