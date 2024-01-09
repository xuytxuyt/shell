using ApproxOperator

include("import_Scordelis_Lo_roof.jl")
ndiv = 1
elements, nodes = import_roof_mix("msh/scordelislo_"*string(ndiv)*".msh");