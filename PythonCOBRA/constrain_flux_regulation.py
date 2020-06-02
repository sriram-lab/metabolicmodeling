"""
constrain_flux_regulation integrates gene expression data into the metabolic model and outputs a flux distribution
that is constrained by gene expression data.
"""

import cobra
import gurobipy as gp
import numpy as np
import scipy.sparse as sp

def constrain_flux_regulation(modelFileName, geneExpressionObj):
    """

    :param modelFileName:
    :param geneExpressionObj:
    :return:
    """

    # Read in common file types for COBRA models
    if modelFileName.split('.')[-1] is ('sbml' or 'xml'):
        cobra_model = cobra.io.read_sbml_model(modelFileName)
    elif modelFileName.split('.')[-1] is 'mat':
        cobra_model = cobra.io.load_matlab_model(modelFileName)

    # Create a model structure
    model = gp.Model("ConstrainedModel")

    # Build the A matrix from the stoichiometric matrix
    S = cobra.util.create_stoichiometric_matrix(cobra_model)
    A = sp.csr_matrix(S, (len(cobra_model.metabolites), len(cobra_model.reactions)))
    lb =
    ub = 
    obj = 
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
    

    
