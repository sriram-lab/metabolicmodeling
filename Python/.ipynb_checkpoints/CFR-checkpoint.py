"""
constrain_flux_regulation integrates gene expression data into the metabolic model and outputs a flux distribution
that is constrained by gene expression data.

@author: Scott Campit
"""

import cobra
import gurobipy as gp
from gurobipy import GRB
import numpy as np
import scipy.sparse as sp

def makeGurobiModel(modelFileName):
    """
    makeGurobiModel creates a Gurobi object from a COBRA model
    """
    # Read in common file types for COBRA models
    if modelFileName.split('.')[-1] == ('sbml' or 'xml'):
        cobra_model = cobra.io.read_sbml_model(modelFileName)
    elif modelFileName.split('.')[-1] == 'mat':
        cobra_model = cobra.io.load_matlab_model(modelFileName)
    print(cobra_model)
    # Begin constructing Gurobi model
    model = gp.Model("COBRA")
    #reactions = model.addVars(cobra_model.reactions, name='Reactions')
    #objective_coefficients = model.addVars()
    #model.setObjective(reactions['biomass_objective'], GRB.MAXMIMIZE)

def constrain_flux_regulation(modelFileName, geneExpressionObj):
    """

    :param modelFileName:
    :param geneExpressionObj:
    :return:
    """

    # Create a model structure
    model = gp.Model()

    # Build the A matrix from the stoichiometric matrix
    S = cobra.util.create_stoichiometric_matrix(cobra_model)
    A = sp.csr_matrix(S, (len(cobra_model.metabolites), len(cobra_model.reactions)))
    lb = None
    ub = None
    obj = None
    # Build the rhs vector
    rhs = np.zeros(np.shape(cobra_model.reactions))

    # Specify additional model parameters
    model.Params.OutputFlag = 0
    model.Params.FeasibilityTol = 1E-9
    model.Params.OptimalityTol = 1E-9
    model.addVar(lb=model.get('lb'),
              ub=model.get('ub'),
              obj=model.get('c'),
              rhs=model.get('b'))
    
if __name__ == "__main__":
    path = "/home/scampit/Data/Reconstructions/RECON1/recon1.mat"
    makeGurobiModel(path)
    
