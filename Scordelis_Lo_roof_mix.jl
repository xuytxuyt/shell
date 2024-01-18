using ApproxOperator

include("import_Scordelis_Lo_roof.jl")
ndiv = 11
elements, nodes = import_roof_mix("msh/scordelislo_"*string(ndiv)*".msh",ndiv-1);