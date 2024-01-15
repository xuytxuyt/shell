
import BenchmarkExample: BenchmarkExample

n = 101
BenchmarkExample.ScordelisLoRoof.generateMsh("./msh/scordelislo_"*string(n)*".msh", transfinite = n)