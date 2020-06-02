"""
Perturbation analyses
This script performs several perturbation analyses using COBRA-based metabolic models.

@author: Scott Campit
"""

import cobra

def geneKnockOut(path2model):
    """
    geneKnockOut returns the objective values resulting from single gene knockouts in a COBRA model.

    :param path2model: A string denoting the path to a metabolic reconstruction.
    :return geneKO: A list of single gene knockout objective values.
    """
    if modelFileName.split('.')[-1] is ('sbml' or 'xml'):
        cobra_model = cobra.io.read_sbml_model(path2model)
    elif modelFileName.split('.')[-1] is 'mat':
        cobra_model = cobra.io.load_matlab_model(path2model)
    elif modelFileName.split('.')[-1] is 'json':
        cobra_model = cobra.io.load_json_model(path2model)

    geneKO = []
    for gene in cobra_model.genes:
        tmp = cobra_model.copy()
        tmp2 = tmp.genes.get_by_id(gene.id).knockout
        geneKO.append(tmp2.optimize().objective_value)
    return geneKO

def reactionKnockOut(path2model):
    """
    reactionKnockOut returns the objective values resulting from single reaction knockouts in a COBRA model.

    :param path2model: A string denoting the path to a metabolic reconstruction.
    :return rxnKO: A list of single reaction knockout objective values.
    """
    if modelFileName.split('.')[-1] is ('sbml' or 'xml'):
        cobra_model = cobra.io.read_sbml_model(path2model)
    elif modelFileName.split('.')[-1] is 'mat':
        cobra_model = cobra.io.load_matlab_model(path2model)
    elif modelFileName.split('.')[-1] is 'json':
        cobra_model = cobra.io.load_json_model(path2model)

    rxnKO = []
    for rxn in cobra_model.rxns:
        tmp = cobra_model.copy()
        tmp2 = tmp.reactions.rxn.bounds = 0, 0
        rxnKO.append(tmp2.optimize().objective_value)
    return rxnKO

def mediumPerturbationAnalyses(path2models, mediumComponents):
    """

    :param path2models:
    :param mediumComponents:
    :return:
    """
    pass