"""
reconstruction.py performs operations that are useful for Genome-Scale Metabolic Reconstructions.

@author: Scott Campit
"""

import pandas as pd

def createGPR(df, id='symbol'):
    """
    createGPR creates a string representation that can be used as a map for the Gene-Protein-Reaction relationship.

    :param df: A pandas dataframe of gene identifiers
    :return:   A pandas dataframe of gene identifiers + the GPR associations (assuming OR logic only)
    """

    gpr = list()
    if id is 'symbol':
        df['AllSymbols'] = df[['Gene symbol', 'Gene symbol synonyms']].astype(str).agg('; '.join, axis=1)
    for row in df['AllSymbols']:
        geneList = row.split("; ")
        #print(geneList)
        gpr.append(' or '.join(['(' + x + ')' for x in geneList]))
    df['GPR'] = gpr
    return df


if __name__ == "__main__":
    fileName='~/Data/Mappings/MetabolicModelMaps/metabolic_map.xlsx'
    sheet='Genes'
    df = pd.read_excel(fileName, sheet_name=sheet)
    gpr_df = createGPR(df)
    gpr_df.to_csv('gpr.csv')
    
        

        
