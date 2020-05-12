# Inferring metabolic fluxes from time-course metabolomics data using Dynamic Flux Analysis (DFA) 
Metabolism is highly dynamic, and is tied to essential cellular processes such as cellular differentiation and proliferation. However, it is challenging to measure the activity of every reaction in a cell. Time-course metabolomics is a useful measure for tracing metabolite dynamics, but most approaches that utilize this data set examines a small part of the metabolic network. 

**Dynamic flux analysis (DFA)** addresses this issue by fitting time-course metabolomics with genome-scale metabolic models to infer metabolic activity throughout the entire metabolic network. 

*Advantages of DFA include*:
  * **Flexibility of data type:** DFA can use both absolute and relative metabolomics values. 
  * **Handling large and small data:** DFA can use as little as two time points, while other models require several time points to fit statistical or genome-scale model. 
  * **Ease of parameter tuning:** DFA models are easily tuned by a single parameter **kappa**, which is a weight proportional to the optimization strength of the objective function.

## Installation
DFA was implemented in [MATLAB](https://www.mathworks.com/products/matlab.html) (recommended version: 2018+), and requires the [Gurobi Mathematical Programming Solver](https://www.gurobi.com/) (recommended version: 2018+). The file `dfa.m` or the livescript `dfa.mlx` in this repository are used to run DFA. 

## Usage
The function `dfa` has the following syntax:
```MATLAB
[dfaModel, solution] = DFA(model, metabolomics, params)
```

**INPUTS:**
  * `model`:             A `structure` of the genome-scale metabolic model in COBRA format. 
  * `metabolomics`:      A `structure` of the metabolomics data. It requires the following `fields`:
    * `data`:            A numerical array of the time course metabolomics data. Rows correspond to metabolites, columns correspond to individual time points. 
    * `position`:        A numerical array mapping the metabolite to compartments in the metabolic model. Rows correspond to metabolites, columns correspond to individual compartments.
  * `params (optional)`: A `structure` containing parameters and hyperparameters. It takes the following `fields`:
    * `kappa`:           A scalar denoting the weight to penalize the metabolomics data. The default value is 1.
    * `kappa2`:          A scalar denoting the weights to minimize the sum of fluxes using parsimonious Flux Balance Analysis (pFBA). The default value is 1E-6.
    * `norm`:            A string denoting which function to use for normalizing the flux activity coefficients. There are three possible arguments. The default argument is `None`:
      * `None`:          No normalization for the flux activity coefficients
      * `MAV`:           Normalizing the flux activity coefficients from the maximum absolute value
      * `Quantile`:      Quantile normalization of the flux activity coefficients
    * `pfba`:            A Boolean that determines whether to run pFBA on or off.

**OUTPUTS:**
  * `dfaModel`:          A `structure` of the genome-scale metabolic model in COBRA format fitted with the flux activity coefficients computed from the time-course metabolomics.
  * `solution`:          A `structure` containing solutions from the Gurobi solver, including the metabolic fluxes.

## Contributing
Contributions are welcome! Please read the contributions guide to get started. Also feel free to submit bugs, feature requests, and pull requests.

Additionally, you can support development for DFA by citing the original manuscript: 

Chandrasekaran, Sriram, Jin Zhang, Zhen Sun, Li Zhang, Christian A. Ross, Yu-Chung Huang, John M. Asara, Hu Li, George Q. Daley, and James J. Collins. ["Comprehensive mapping of pluripotent stem cell metabolism using dynamic genome-scale network modeling."](https://www.cell.com/cell-reports/fulltext/S2211-1247(17)31027-6) Cell reports 21, no. 10 (2017): 2965-2977.

## Publications using Dynamic Flux Analysis
**Stem cell metabolism**
1. Chandrasekaran, Sriram, Jin Zhang, Zhen Sun, Li Zhang, Christian A. Ross, Yu-Chung Huang, John M. Asara, Hu Li, George Q. Daley, and James J. Collins. ["Comprehensive mapping of pluripotent stem cell metabolism using dynamic genome-scale network modeling."](https://www.cell.com/cell-reports/fulltext/S2211-1247(17)31027-6) Cell reports 21, no. 10 (2017): 2965-2977.
2. Shen, Fangzhou, Camden Cheek, and Sriram Chandrasekaran. ["Dynamic network modeling of stem cell metabolism."](https://link.springer.com/protocol/10.1007%2F978-1-4939-9224-9_14) In Computational stem cell biology, pp. 305-320. Humana, New York, NY, 2019.

**Cancer metabolism**
1. Nelson, Barbara S., Lin Lin, Daniel M. Kremer, Cristovão M. Sousa, Cecilia Cotta-Ramusino, Amy Myers, Johanna Ramos et al. ["Tissue of origin dictates GOT1 dependence and confers synthetic lethality to radiotherapy."](https://cancerandmetabolism.biomedcentral.com/articles/10.1186/s40170-019-0202-2) Cancer & Metabolism 8, no. 1 (2020): 1-16.
2. Campit, Scott, and Sriram Chandrasekaran. ["Inferring Metabolic Flux from Time-Course Metabolomics."](https://link.springer.com/protocol/10.1007%2F978-1-0716-0159-4_13) In Metabolic Flux Analysis in Eukaryotic Cells, pp. 299-313. Humana, New York, NY, 2020.
