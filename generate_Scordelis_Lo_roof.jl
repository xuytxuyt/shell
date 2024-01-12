
import BenchmarkExample: BenchmarkExample

n = 51
BenchmarkExample.ScordelisLoRoof.generateMsh("./msh/scordelislo_"*string(n)*".msh", transfinite = n)