
using ApproxOperator, GLMakie

import BenchmarkExample: BenchmarkExample
import Gmsh: gmsh

n = 56
# lc = 0.1
# lc = 0.3
# filename = "patchtest"
# BenchmarkExample.ScordelisLoRoof.generateMsh("./msh/scordelislo_"*string(n)*".msh", transfinite = n)
BenchmarkExample.SphericalShell.generateMsh("./msh/sphericalshell_"*string(n)*".msh", transfinite = n)
# BenchmarkExample.PatchTestThinShell.generateMsh("./msh/"*filename*".msh",lc = lc)
