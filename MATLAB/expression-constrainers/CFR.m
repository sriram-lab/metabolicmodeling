function [ConstrainedModel, solution] = CFR(model, hyperparams, onreactions, offreactions)
%% CFR Using transcriptomics to compute context specific metabolic fluxes
% The *Constrained Flux Regulation* (|CFR|) function is a hybrid of the original 
% iMAT algorithm (<https://academic.oup.com/bioinformatics/article/26/24/3140/290045 
% Zur et al., 2010>) and the PROM algorithm (<https://www.pnas.org/content/107/41/17845.long 
% Chandrasekaran et al., 2010>). An informal description of the algorithm is provided 
% below:
%% 
% * Reactions associated with upregulated genes have their flux values maximized. 
% * Reactions associated with downregulated genes have their flux values minimized. 
% * If the |pfba| flag in the |hyperparams| structure is set to |true|, the 
% sum of fluxes through the entire metabolic network is minimized. This provides 
% a unique flux distribution.
%% 
% A more formal description of the algorithm can be found in <https://drive.google.com/file/d/1XiBtrznWdwlXj9dHRqrFsKKZ0O9sPdFZ/view 
% this book chapter>.
% 
% *INPUTS*
% 
% |model|: A structure of the genome-scale metabolic model in COBRA format.
% 
% |hyperparams|: A structure containing parameters for constraining the flux 
% solution. You can read <https://drive.google.com/file/d/1XiBtrznWdwlXj9dHRqrFsKKZ0O9sPdFZ/view 
% this book chapter> for more details about individual parameters.
% 
% |onreactions|: A cell array containing gene symbols or BiGG reaction names 
% that are upregulated.
% 
% |offreactions|: A cell array containing gene symbols or BiGG reaction names 
% that are downregulated.
% 
% *OUTPUTS*
% 
% |ConstrainedModel:| A structure of the transcriptomics-constrained genome-scale 
% metabolic model.
% 
% |solution:| A Gurobi structure containing the linear programming solution.
    
    % Unpack params
    epsilon         = hyperparams.eps;
    kappa           = hyperparams.kap;
    rho             = hyperparams.rho;
    mode            = hyperparams.mode;
    epsilon2        = hyperparams.eps2;
    minfluxflag     = hyperparams.pfba;
    kappa2          = hyperparams.kap2;
    
    % Set upregulated / downregulated to empty cell array if empty
    if (~exist('onreactions')) || (isempty(onreactions))
        onreactions = {};
    end
    if (~exist('offreactions')) || (isempty(offreactions))
        offreactions = {};
    end
    
    % Find reactions associated with gene symbols
    if (~exist('mode','var')) || (isempty(mode))
        mode = 1;
    elseif mode == 1 
        [~, ~,  onreactions, ~] =  deleteModelGenes(model, cellstr(onreactions));
        [~, ~, offreactions, ~] =  deleteModelGenes(model, cellstr(offreactions));
    end
    onreactions = unique(onreactions);
    offreactions = unique(offreactions);
    
    % Set penalty for downregulated reaction
    if (~exist('kappa','var')) || (isempty(kappa))
        kappa = ones(size(offreactions));
    end
    if numel(kappa) == 1
        kappa = ones(size(offreactions));
    end
    
    % Set weight for upregulated reaction
    if (~exist('rho','var')) || (isempty(rho))
        rho = ones(size(onreactions));
    end
    if numel(rho) == 1
        rho = ones(size(onreactions));
    end
    
    % Set the minimum flux for upregulated reactions
    if (~exist('epsilon','var')) || (isempty(epsilon))
        epsilon = ones(size(onreactions)) .* 1E-3;
    end
    if numel(epsilon) == 1
        epsilon = repmat(epsilon, size(onreactions));
    end
    
    % Set the right hand side for downregulated reactions
    if (~exist('epsilon2','var')) || (isempty(epsilon2))
        epsilon2 = zeros(size(offreactions));
    end
    
    % Set Parsimonious Flux Balance Analysis 
    if isempty(minfluxflag)
        minfluxflag = true;
    end
    if minfluxflag == true
        kappa    = [kappa(:); ones(size(setdiff(model.rxns, offreactions))) * 1E-4];
        epsilon2 = [epsilon2(:); zeros(size(setdiff(model.rxns, offreactions)))];
        offreactions = [offreactions(:); setdiff(model.rxns, offreactions)];
    end
    
    % Setup Gurobi model
    tmp     = model;
    tmp.A   = tmp.S;
    tmp.obj = tmp.c;
    tmp.rhs = tmp.b;
    if (exist('model.csense','var')) && (~isempty(model.csense))
        tmp.sense = model.csense;
        tmp.sense(ismember(model.sense,'E')) = '=';
        tmp.sense(ismember(model.sense,'L')) = '<';
        tmp.sense(ismember(model.sense,'G')) = '>';
    else
        tmp.sense = repmat('=', [size(model.S, 1), 1]);
    end
    
    tmp.lb         = model.lb;
    tmp.ub         = model.ub;
    tmp.vtype      = repmat('C', size(model.S, 2), 1);
    tmp.modelsense = 'max';
    M              = 10000;
    
    % RESET OBJECTIVE TO BIOMASS REACTION
    %tmp.c          = zeros(size(tmp.c));
    %biomassPos     = contains(tmp.rxns, 'biomass');
    %tmp.c          = biomassPos;
    
    % Maximize the flux through reactions that have upregulated genes
    function tmp1 = max_flux_through_oneactions(onreactions, tmp, epsilon, rho, M)
        tmp1 = tmp;
        for i = 1:length(onreactions)
            onpos = find(ismember(tmp1.rxns, onreactions(i)));
            
            % Set reaction in forward reaction to be active
            newMet = size(tmp1.A, 1) + 1; newRxn = size(tmp1.A, 2) + 1;
            
            tmp1.A(newMet, onpos)   = 1;
            tmp1.A(newMet, newRxn)  = -(1 * epsilon(i) + M);
            tmp1.rhs(newMet)        = -M;
            tmp1.sense(newMet)      = '>';
            tmp1.vtype(newRxn)      = 'B';
            tmp1.obj(newRxn)        = 1 * rho(i);
            tmp1.lb(newRxn)         = 0; 
            tmp1.ub(newRxn)         = 1;
            
            % Set reaction in backward reaction to be active
            newMet = size(tmp1.A, 1) + 1; newRxn = size(tmp1.A, 2) + 1;
            
            tmp1.A(newMet, onpos)   = 1;
            tmp1.A(newMet, newRxn)  = (1*epsilon(i) + M);
            tmp1.rhs(newMet)        = M;
            tmp1.sense(newMet)      = '<';
            tmp1.vtype(newRxn)      = 'B';
            tmp1.obj(newRxn)        = 1 * rho(i);
            tmp1.lb(newRxn)         = 0; 
            tmp1.ub(newRxn)         = 1;
        end
    end
    
    function tmp2 = min_flux_through_offeactions(offreactions, tmp, epsilon2, kappa)
        
        tmp2 = tmp;
        % Minimize flux through reactions that have downregulated genes
        for j = 1:length(offreactions)
            offpos = find(ismember(tmp.rxns, offreactions(j)));
            
            % xi + si >= -eps2
            %      si >= 0
            % rho(ri + si)
            newMet = size(tmp2.A, 1) + 1;
            newRxn = size(tmp2.A, 2) + 1;
            
            tmp2.A(newMet, offpos) = 1;
            tmp2.A(newMet, newRxn) = 1;
            tmp2.rhs(newMet)       = -1 * epsilon2(j);
            tmp2.sense(newMet)     = '>';
            tmp2.vtype(newRxn)     = 'C';
            tmp2.lb(newRxn)         = 0; 
            tmp2.ub(newRxn)         = 1000;
            tmp2.obj(newRxn)       = -1 * kappa(j);
            
            % xi - ri <= eps2
            %      ri >= 0
            newMet = size(tmp2.A, 1) + 1;
            newRxn = size(tmp2.A, 2) + 1;
            
            tmp2.A(newMet, offpos) = 1;
            tmp2.A(newMet, newRxn) = -1;
            tmp2.rhs(newMet)       = epsilon2(j);
            tmp2.sense(newMet)     = '<';
            tmp2.vtype(newRxn)     = 'C';
            tmp2.lb(newRxn)         = 0; 
            tmp2.ub(newRxn)         = 1000;
            tmp2.obj(newRxn)       = -1 * kappa(j);
        end
    end
    
    % Set params for Gurobi
    params.outputflag          = 0;
    params.Threads             = 2;
    params.Seed                = 314;
    params.NumericFocus        = 3;
    
    up  = max_flux_through_oneactions(onreactions, tmp, epsilon, rho, M);
    up_and_down = min_flux_through_offeactions(offreactions, up, epsilon2, kappa); 
    
    % Solve using the Gurobi solver
    ConstrainedModel           = up_and_down;
    solution                   = gurobi(ConstrainedModel, params);
    grate                      = solution.x(logical(ConstrainedModel.c));
    
    % Check and see if a biomass constraint was added onto the model. If
    % too constrained, square it until a feasible flux solution is
    % obtained.
    %if model.ub(logical(ConstrainedModel.c)) < 1000
    %    while grate == 0
    %        tmp.ub(logical(ConstrainedModel.c)) = tmp.ub(logical(ConstrainedModel.c))^2;
    %        up                 = max_flux_through_oneactions(onreactions,   tmp, epsilon, rho, M);
    %        up_and_down        = min_flux_through_offeactions(offreactions, up, epsilon2, kappa); 
    %        ConstrainedModel   = up_and_down;
    %        solution           = gurobi(ConstrainedModel, params);
    %        grate              = solution.x(logical(ConstrainedModel.c));
    %    end
    %end
    
    k = 1;
    while grate == 0
        epsilon2 = epsilon2 .* 0.1;
        kappa    = kappa    .* 0.1;
        
        up                 = max_flux_through_oneactions(onreactions,   tmp, epsilon, rho, M);
        up_and_down        = min_flux_through_offeactions(offreactions, up, epsilon2, kappa); 
        ConstrainedModel   = up_and_down;
        solution           = gurobi(ConstrainedModel, params);
        grate              = solution.x(logical(ConstrainedModel.c));
        
        if k > 100
            break
        else
            k = k + 1;
        end   
    end
end