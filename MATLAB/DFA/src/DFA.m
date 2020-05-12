function [dynamicModel, solution] = DFA(model, metabolomics, params)
%% DFA Dynamic flux analysis for time-course metabolomics data
% |dfa| computes metabolic fluxes using time-course metabolomics data.
% 
% *INPUTS*
%% 
% # |model|:                A geneome-scale metabolic model in COBRA format
% # |file|:                  An excel spreadsheet containing the metabolomics 
% dataset. 
% # |sheetname|:       The sheet name in the Excel spreadsheet. The default 
% sheet is 'Sheet1'
% # |kappa|:                The relative weight for the consistency of the flux 
% distribution compared to the weight for maximizing the objective function. The 
% default value is 1
% # |kappa2|:              The relative weight for minimizing the sum of absolute 
% flux values of other reactions. The default value is 1E-3
%% 
% *OUTPUTS*
%% 
% # |dynamicModel|:  A metabolic model with pseudoreactions that correspond 
% to the time-course metabolomics data
% # |solution|:          The growth rate objective value from single gene KO
    
    % Default parameters
    if (~exist('params.kappa', 'var')) || (isempty(params.kappa))
        params.kappa = 1;
    end
    if (~exist('params.kappa2', 'var')) || (isempty(params.kappa2))
        params.kappa2 = 1E-3;
    end
    if (~exist('params.norm', 'var')) || (isempty(params.norm))
        params.norm = 'None';
    end
    
    % Initialize Gurobi model parameters
    tmpModel               = model;
    tmpModel.A             = model.S;
    tmpModel.obj           = model.c;
    tmpModel.rhs           = model.b;
    tmpModel.sense         = repmat('=', [size(model.S, 1), 1]);
    tmpModel.lb            = model.lb;
    tmpModel.ub            = model.ub;
    tmpModel.vtype         = repmat('C', [size(model.S, 2), 1]);
    tmpModel.modelsense    = 'max';
    tmpModel.original_size = length(tmpModel.rxns);
    
    % Compute the linear flux activity coefficient
    slope = zeros(size(metabolomics.data, 1), 1);
    for m = 1:size(metabolomics.data, 1)
        weight   = polyfit(1:size(metabolomics.data, 2), metabolomics.data(m, :), 1);
        slope(m) = weight(1) / weight(2);
    end
    
    % Normalize the flux activity coefficients
    switch params.norm
        case 'None'
            slope                      = slope;
        case 'MAV'
            slope(:, 1)                = slope(:, 1) ./ max(abs(slope(:, 1)));
        case 'Quantile'
            slope(:, 1)                = slope(:, 1) ./ quantilenorm(abs(slope(:, 1)), 'MEDIAN', true);
            slope(slope(:, 1) > 1, 1)  = log10(slope(slope(:, 1) > 1, 1)) + 1;
            slope(slope(:, 1) < -1, 1) = -log10(-1 * slope(slope(:, 1) < -1, 1)) - 1;
    end
    
    % Create the right hand side and constrain the model using the flux
    % activity coefficients
    constrainedModel = tmpModel;
    epsilon          = zeros(size(slope(:, 1)));
    for m = 1:size(metabolomics.positions, 1)
        metabolite = (metabolomics.positions ~= 0);
        if ~isempty(metabolite)
            positions                  = metabolomics.positions(m, :);
            positions(positions == 0)  = [];
            epsilon(m)                 = slope(m); 
    
            alpha = size(constrainedModel.A, 2) + 1;
            beta  = alpha + 1;
            gamma = size(constrainedModel.A, 1) + 1;
            constrainedModel.rhs(positions) = epsilon(m);
            constrainedModel.rhs(isnan(constrainedModel.rhs)) = 0;
    
            % metabolite constraint = tolerance - alpha + beta
            constrainedModel.A(positions, alpha) = 1;
            constrainedModel.A(positions, beta)  = -1;
            constrainedModel.vtype(alpha)        = 'C';
            constrainedModel.lb(beta)            = 0;
            constrainedModel.lb(alpha)           = 0;
            constrainedModel.ub(beta)            = 1000;
            constrainedModel.ub(alpha)           = 1000;
            constrainedModel.obj(beta)           = -1 * params.kappa; 
            constrainedModel.obj(alpha)          = -1 * params.kappa;
            constrainedModel.c(alpha)            = -1 * params.kappa;
            constrainedModel.c(beta)             = -1 * params.kappa;
            constrainedModel.rxns{alpha}         = 'alpha';
            constrainedModel.rxns{beta}          = 'beta';
            constrainedModel.mets{gamma}         = 'pmet';
            constrainedModel.sense(alpha)        = '=';
            constrainedModel.sense(beta)         = '=';
            
        end
    end
    
    % Minimize the sum of absolute flux values for all other reactions
    % using kappa2 (pFBA)
    for r = 1:length(constrainedModel.rxns)
        if constrainedModel.c(r) == 0
            
            % xi + si >= -eps2
            %      si >= 0
            % rho(ri + si)
            pmet                           = size(constrainedModel.A, 1) + 1;
            prxn                           = size(constrainedModel.A, 2) + 1;
            constrainedModel.A(pmet, r)    = 1;
            constrainedModel.A(pmet, prxn) = 1;
            constrainedModel.rhs(pmet)     = -0;
            constrainedModel.sense(pmet)   = '>';
            constrainedModel.lb(prxn)      = 0;
            constrainedModel.ub(prxn)      = 1000;
            constrainedModel.obj(prxn)     = -1 * params.kappa2;
            constrainedModel.c(prxn)       = -1 * params.kappa2;
            constrainedModel.rxns{prxn}    = 'pseudorxn';
            % xi - ri <= eps2
            %      ri >= 0
            pmet                           = size(constrainedModel.A, 1) + 1;
            prxn                           = size(constrainedModel.A, 2) + 1;
            constrainedModel.A(pmet, r)    = 1;
            constrainedModel.A(pmet, prxn) = -1;
            constrainedModel.rhs(pmet)     = 0;
            constrainedModel.sense(pmet)   = '<';
            constrainedModel.lb(prxn)      = 0;
            constrainedModel.ub(prxn)      = 1000;
            constrainedModel.obj(prxn)     = -1 * params.kappa2;
            constrainedModel.rxns{prxn}    = 'pseudorxn';
        end
    end
    
    constrainedModel.vtype     = repmat('C', size(constrainedModel.A, 2), 1);
    dynamicModel               = constrainedModel;
    
    % Set params for Gurobi
    params2.outputflag          = 0;
    params2.Threads             = 2;
    params2.Seed                = 314;
    params2.NumericFocus        = 3;
    solution                    = gurobi(dynamicModel, params2);
    
    try
        solution.v              = solution.x(1:length(model.rxns));
        solution.solverObj      = solution.objval;
        solution.rhs            = slope;
    catch ME
        warning('Infeasible solution obtained for a metabolic model')
    end
end