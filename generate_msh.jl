
import BenchmarkExample: BenchmarkExample

n = 1
# BenchmarkExample.ScordelisLoRoof.generateMsh("./msh/scordelislo_"*string(n)*".msh", transfinite = n)
BenchmarkExample.SphericalShell.generateMsh("./msh/sphericalshell_"*string(n)*".msh", transfinite = n)