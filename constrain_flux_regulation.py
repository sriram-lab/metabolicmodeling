"""

"""

import cobra
import gurobi

def constrain_flux_regulation(model, geneExpressionObj):
    lp = gurobipy.Model()
    lp.Params.OutputFlag = 0
    lp.Params.FeasibilityTol = 1E-9
    lp.Params.OptimalityTol = 1E-9
    lp.addVar(lb=model.get('lb'),
              ub=model.get('ub'),
              obj=model.get('c'),
              rhs=model.get('b'))
    

    
