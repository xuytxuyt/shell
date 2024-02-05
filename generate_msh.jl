
import BenchmarkExample: BenchmarkExample

n = 1
lc = 0.1
# BenchmarkExample.ScordelisLoRoof.generateMsh("./msh/scordelislo_"*string(n)*".msh", transfinite = n)
# BenchmarkExample.SphericalShell.generateMsh("./msh/sphericalshell_"*string(n)*".msh", transfinite = n)
BenchmarkExample.PatchTestThinShell.generateMsh("./msh/patchtest.msh",lc = lc)