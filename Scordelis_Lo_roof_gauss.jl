using ApproxOperator

include("import_Scordelis_Lo_roof.jl")
ndiv = 1
elements, nodes = import_roof_gauss("msh/scordelislo_"*string(ndiv)*".msh");