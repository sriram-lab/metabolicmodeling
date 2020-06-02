function medium_model = addMediumConstraints(model, queried_medium, path_to_file)
%% ADDMEDIUMCONSTRAINTS Add medium constraints to metabolic model
% |addMediumConstraints| sets the lower bound of the metabolic model to reflect 
% different medium conditions. 
% 
% *Note:* The medium condition file is assumed to be an Excel sheet, where each 
% sheet name is the name of the medium component. The following details should 
% be included in your medium map:
%% 
% # The medium component common name
% # The default lowerbound you intend to set for the metabolic model (Default 
% lowerbound = -1.0)
% # The scaling factor (alpha), which I have defined as the proportion of the 
% new medium component with respect to the default medium component that the metabolic 
% model's lower bounds are set to. (Default base medium for RECON1 = RPMI)
% # The adjusted lowerbound (default lb x alpha)
% # The BiGG Reaction ID
% # The BiGG Metabolite ID
%% 
% *INPUTS*
%% 
% * |model|: A MATLAB structure that contains the lower bounds for the constraint-based 
% model.
% * |path_to_file|: A string denoting the path to the medium sheet used to set 
% the lower bounds. 
% * |queried_medium|: A string denoting the medium condition being simulated.
%% 
% *OUTPUT*
%% 
% * |medium_model|: A MATLAB structure that contains the constrained lower bounds.
    
    if nargin < 3
        path_to_file = 'FINAL_MEDIUM_MAP.xlsx';
    end
    
    medium_map = path_to_file;
    medium_model = model;
    
    % Support for MATLAB versions less than 2018
    if verLessThan('matlab', '9.6.0.1072779')
        [~, sheetNames] = xlsfinfo(medium_map);
        for sheets = 1:length(sheetNames)
            
            % If the sheet name matches the queried medium, set the lower
            % bounds using the BiGG Reaction IDs
            if ismember(string(sheetNames(sheets)), queried_medium)
                [adjustedLB, rxn_ids] = xlsread(medium_map, string(sheetNames{sheets}));
                rxn_ids(1, :) = [];
                rxn_ids(:, 1) = [];
                adj_lb = adjustedLB(:, 2);
                
                for rxn=1:length(rxn_ids)
                    pos = find(ismember(string(medium_model.rxns), string(rxn_ids(rxn, 5))));
                    medium_model.lb(pos) = adj_lb(rxn);
                end
            
            % If the sheet name cannot be found, set the lower bounds to
            % RPMI conditions.
            elseif ismember({'nan'}, queried_medium)
                [adjustedLB, rxn_ids] = xlsread(medium_map, 'RPMI');
                rxn_ids(1, :) = [];
                rxn_ids(:, 1) = [];
                adj_lb = adjustedLB(:, 2);
                
                for rxn=1:length(rxn_ids)
                    pos = find(ismember(string(medium_model.rxns), string(rxn_ids(rxn, 5))));
                    medium_model.lb(pos) = adj_lb(rxn);
                end
            end
        end 
    
    % Support for newer versions of MATLAB
    else
        [~, sheetNames] = xlsfinfo(medium_map);
        for sheets = 1:length(sheetNames)
            
            % If the sheet name matches the queried medium, set the lower
            % bounds using the BiGG Reaction IDs
            if ismember(string(sheetNames(sheets)), string(queried_medium))
                dataArray = readcell(medium_map, 'Sheet', string(sheetNames(sheets)));
                dataArray(1,:) = [];
                dataArray(:,1) = [];
                adj_lb = cell2mat(dataArray(:, 2));
                
                for rxn = 1:size(dataArray, 1)
                    pos = find(ismember(string(medium_model.rxns), string(dataArray(rxn, 4))));
                    medium_model.lb(pos) = adj_lb(rxn);
                end
            
            % If the sheet name cannot be found, set the lower bounds to
            % RPMI conditions.
            elseif ismember({'nan'}, queried_medium)
                dataArray = readcell(medium_map, 'Sheet', 'RPMI');
                dataArray(1,:) = [];
                dataArray(:,1) = [];
                adj_lb = cell2mat(ataArray(:, 2));
                
                for rxn = 1:size(dataArray, 1)
                    pos = find(ismember(string(medium_model.rxns), string(dataArray(rxn, 4))));
                    medium_model.lb(pos) = adj_lb(rxn);
                end    
            end
        end
    end   
end