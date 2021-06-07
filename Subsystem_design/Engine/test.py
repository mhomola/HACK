import cantera
gas = cantera.Solution('gri30.yaml')
h2o = cantera.PureFluid('liquidvapor.yaml', 'water')