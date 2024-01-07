
import BenchmarkExample: BenchmarkExample

n = 10
BenchmarkExample.ScordelisLoRoof.generateMsh("./msh/scordelislo_"*string(n)*".msh", transfinite = n)