function [ConstrainedModel, solution] = CFR2(model, hyperparams, reactions)
%% CFR2 Constrain metabolic fluxes using the modified iMAT algorithm and PROM lowerbound constraint
% 
    
    % Unpack params
    epsilon         = hyperparams.epsilon;
    kappa           = hyperparams.kappa;
    rho             = hyperparams.rho;
    mode            = hyperparams.mode;
    epsilon2        = hyperparams.epsilon2;
    minfluxflag     = hyperparams.pfba;
    
    
    % Set the right hand side for downregulated reactions
    if (~exist('epsilon2','var')) || (isempty(epsilon2))
        epsilon2 = zeros(size(reactions));
    end
    
    % Set Parsimonious Flux Balance Analysis 
    if isempty(minfluxflag)
        minfluxflag = true;
    end
    if numel(minfluxflag) == 1
        kappa    = [kappa(:); ...
            ones(size(setdiff(model.rxns, reactions))) * 1E-6];
        epsilon2 = [epsilon2(:); ...
            zeros(size(setdiff(model.rxns, reactions)))];
        offreactions = [reactions(:); ...
            setdiff(model.rxns, reactions)];
    end
    
    % Setup Gurobi model
    tmp     = model;
    tmp.A   = tmp.S;
    tmp.obj = tmp.c;
    tmp.rhs = tmp.b;
    if (exist('tmp.csense','var')) && (~isempty(tmp.csense))
        tmp.sense = tmp.csense;
        tmp.sense(ismember(tmp.sense,'E')) = '=';
        tmp.sense(ismember(tmp.sense,'L')) = '<';
        tmp.sense(ismember(tmp.sense,'G')) = '>';
    else
        tmp.sense = repmat('=', [size(tmp.S, 1), 1]);
    end
    
    tmp.lb         = tmp.lb;
    tmp.ub         = tmp.ub;
    tmp.vtype      = repmat('C', size(tmp.S, 2), 1);
    tmp.modelsense = 'max';
    M              = 10000;
    
    % RESET OBJECTIVE TO BIOMASS REACTION
    %tmp.c          = zeros(size(tmp.c));
    %biomassPos     = contains(tmp.rxns, 'biomass');
    %tmp.c          = biomassPos;
    
    % Maximize the flux through reactions that have upregulated genes
    for i = 1:length(reactions)
        onpos = find(ismember(tmp.rxns, reactions(i)));
        
        % Set reaction in forward reaction to be active
        newMet = size(tmp.A, 1) + 1; newRxn = size(tmp.A, 2) + 1;
        
        tmp.A(newMet, onpos)   = 1;
        tmp.A(newMet, newRxn)  = -(1 * epsilon(i) + M);
        tmp.rhs(newMet)        = -M;
        tmp.sense(newMet)      = '>';
        tmp.vtype(newRxn)      = 'B';
        tmp.obj(newRxn)        = 1 * rho(i);
        tmp.lb(newRxn)         = min(tmp.lb); 
        tmp.ub(newRxn)         = max(tmp.ub);
        
        % Set reaction in backward reaction to be active
        newMet = size(tmp.A, 1) + 1; newRxn = size(tmp.A, 2) + 1;
        
        tmp.A(newMet, onpos)   = 1;
        tmp.A(newMet, newRxn)  = (1*epsilon(i) + M);
        tmp.rhs(newMet)        = M;
        tmp.sense(newMet)      = '<';
        tmp.vtype(newRxn)      = 'B';
        tmp.obj(newRxn)        = 1 * rho(i);
        tmp.lb(newRxn)         = min(tmp.lb); 
        tmp.ub(newRxn)         = max(tmp.ub);
    end
    
    % Minimize flux through reactions that have downregulated genes
    for j = 1:length(reactions)
        offpos = find(ismember(tmp.rxns, ...
            reactions(j)));
        
        % xi + si >= -eps2
        %      si >= 0
        % rho(ri + si)
        newMet = size(tmp.A, 1) + 1;
        newRxn = size(tmp.A, 2) + 1;
        
        tmp.A(newMet, offpos) = 1;
        tmp.A(newMet, newRxn) = 1;
        tmp.rhs(newMet)       = -1 * epsilon2(j);
        tmp.sense(newMet)     = '>';
        tmp.vtype(newRxn)     = 'C';
        tmp.lb(newRxn)        = min(tmp.lb);
        tmp.ub(newRxn)        = max(tmp.ub);
        tmp.obj(newRxn)       = -1 * kappa(j);
        
        % xi - ri <= eps2
        %      ri >= 0
        newMet = size(tmp.A, 1) + 1;
        newRxn = size(tmp.A, 2) + 1;
        
        tmp.A(newMet, offpos) = 1;
        tmp.A(newMet, newRxn) = -1;
        tmp.rhs(newMet)       = epsilon2(j);
        tmp.sense(newMet)     = '<';
        tmp.vtype(newRxn)     = 'C';
        tmp.lb(newRxn)        = min(tmp.lb);
        tmp.ub(newRxn)        = max(tmp.ub);
        tmp.obj(newRxn)       = -1 * kappa(j);
    end
    
    % Set params for Gurobi
    params.outputflag          = 0;
    params.Threads             = 2;
    params.Seed                = 314;
    params.NumericFocus        = 3;
    
    % Solve using the Gurobi solver
    ConstrainedModel           = tmp;
    solution                   = gurobi(ConstrainedModel, params);
    %ConstrainedModel.vbasis    = solution.vbasis;
    %ConstrainedModel.cbasis    = solution.cbasis;
    %ConstrainedModel.slack     = solution.slack;
    %ConstrainedModel.pi        = solution.pi;
    %ConstrainedModel.rc        = solution.rc;
    %solution                   = gurobi(ConstrainedModel, params);
end