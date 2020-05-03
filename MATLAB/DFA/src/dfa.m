%% @authors: Sriram Chandrasekaran, Fangzhou Shen, Scott Campit
function [dynamicModelFlux, dynamicmodel] = dfa(model, metabolomics_sheet, ...
        sheetname, kappa, kappa2, normalizationFlag, timeseries)

if (~exist('kappa', 'var')) || (isempty(kappa))
        kappa = 1;
end
if (~exist('kappa2', 'var')) || (isempty(kappa2))
        kappa2 = 1E-3;
end
if (~exist('sheetname', 'var')) || (isempty(sheetname))
        sheetname = 'Sheet1';
end
if (~exist('genedelflag', 'var')) || (isempty(genedelflag))
        genedelflag = 0;
end
if (~exist('rxndelflag', 'var')) || (isempty(rxndelflag))
        rxndelflag = 0;
end

        function gurobi_model = make_gurobi_model(model)
                gurobi_model = model;
                gurobi_model.A = model.S;
                gurobi_model.obj = model.c;
                gurobi_model.rhs = model.b;
                gurobi_model.sense = repmat('=', [size(model.S, 1), 1]);
                gurobi_model.lb = model.lb;
                gurobi_model.ub = model.ub;
                gurobi_model.vtype = repmat('C', [size(model.S, 2), 1]);
                gurobi_model.modelsense = 'max';
                gurobi_model.original_size = length(gurobi_model.rxns);
        end

        function [time_points, model_positions, final_metabolomics] =...
                        preprocess_metabolomics(metabolomics_sheet, sheetname, timeseries)
                
                [values, ~] = xlsread(metabolomics_sheet, sheetname);
                %if size(values, 2) > 14 % shME1
                switch timeseries
                        case '0 to 5'
                                time_points = values(1, 9:20);
                                model_positions = values(2:end, 1:8);
                                metabolomics = values(2:end, 9:20);
                                final_metabolomics = knnimpute(metabolomics);
                        case '1 to 5'
                                time_points = values(1, 12:20);
                                model_positions = values(2:end, 1:8);
                                metabolomics = values(2:end, 12:end);
                                final_metabolomics = knnimpute(metabolomics);
                        case '0 to 1'
                                time_points = values(1, 9:14);
                                model_positions = values(2:end, 1:8);
                                metabolomics = values(2:end, 9:14);
                                final_metabolomics = knnimpute(metabolomics);
                        case '1 to 3'
                                time_points = values(1, 12:17);
                                model_positions = values(2:end, 1:8);
                                metabolomics = values(2:end, 12:17);
                                final_metabolomics = knnimpute(metabolomics);
                        case '3 to 5'
                                time_points = values(1, 15:20);
                                model_positions = values(2:end, 1:8);
                                metabolomics = values(2:end, 15:20);
                                final_metabolomics = knnimpute(metabolomics);
                end
                %else % shNT
                %    time_points = [1, 1, 1, 5, 5, 5];
                %    model_positions = values(2:end,1:8);
                %    metabolomics = values(2:end, 9:end);
                %    final_metabolomics = knnimpute(metabolomics);
                %end
        end

        function slope = compute_fac(final_metabolomics, time_points)
                slope = zeros(size(final_metabolomics, 1), 2);
                for metabolite = 1:size(final_metabolomics, 1)
                        [poly_coef, ~] = polyfit(time_points, final_metabolomics(metabolite, :), 1);
                        slope(metabolite, :) = poly_coef(1) / poly_coef(2);
                end
        end

        function normedSlope = normalize_fac(slope, normalizationFlag)
                switch normalizationFlag
                        case 'maximumAbsoluteValue'
                                normedSlope(:,1) = slope(:,1)./ max(abs(slope(:,1)));
                        case 'Quantile'
                                normedSlope(:, 1) = slope(:,1)./ quantilenorm(abs(slope(:,1)), 'MEDIAN', true);
                                normedSlope(normedSlope(:, 1) > 1, 1) = log10(normedSlope(normedSlope(:, 1) > 1, 1)) + 1;
                                normedSlope(normedSlope(:, 1) < -1, 1) = -log10(-1*normedSlope(normedSlope(:, 1) < -1, 1)) - 1;
                end
        end

        function constrain_model = addMetabolomicConstraints(gurobi_model,...
                        final_metabolomics, model_positions, normedSlope, kappa)
                constrain_model = gurobi_model;
                for metabolite = 1:size(final_metabolomics, 1)
                        metabolite_in_compartment = (model_positions(:, 1) ~= 0);
                        flux_activity_coef = zeros(size(normedSlope(:, 1)));
                        
                        if metabolite_in_compartment(metabolite)
                                compartment_positions = model_positions(metabolite, :);
                                compartment_positions(compartment_positions == 0) = '';
                                flux_activity_coef(metabolite) = normedSlope(metabolite, 1);
                                epsilon = flux_activity_coef(metabolite);
                                
                                rxnpos = compartment_positions;
                                alpha = size(constrain_model.A, 2) + 1;
                                beta = alpha + 1;
                                disp(rxnpos)
                                constrain_model.rhs(rxnpos) = epsilon;
                                
                                % metrow = tol - alpha + beta
                                constrain_model.A(rxnpos, alpha) = 1;
                                constrain_model.A(rxnpos, beta) = -1;
                                constrain_model.vtype(alpha) = 'C';
                                
                                % set si to be positive
                                constrain_model.lb(beta) = 0;
                                constrain_model.ub(beta) = 100;
                                constrain_model.lb(alpha) = 0;
                                constrain_model.ub(alpha) = 100;
                                constrain_model.obj(beta) = -1*kappa;
                                constrain_model.obj(alpha) = -1*kappa;
                                constrain_model.c(alpha) = -1*kappa;
                                constrain_model.c(beta) = -1*kappa;
                                constrain_model.rxns{alpha} = 'alpha';
                                constrain_model.rxns{beta} = 'beta';
                        end
                end
        end

        function dynamicmodel = setDefaultConstraints(constrain_model, kappa2)
                dynamicmodel = constrain_model;
                for rxn = 1:length(dynamicmodel.rxns)
                        if dynamicmodel.c(rxn) == 0
                                % xi + si >= -eps2
                                % si >= 0
                                % rho(ri + si)
                                
                                % constraint 1
                                rxnpos = size(dynamicmodel.A, 1) + 1;
                                alpha = size(dynamicmodel.A, 2) + 1;
                                dynamicmodel.A(rxnpos, rxn) = 1;
                                dynamicmodel.A(rxnpos, alpha) = 1;
                                dynamicmodel.rhs(rxnpos) = 0;
                                dynamicmodel.sense(rxnpos) = '>';
                                
                                % set si to be positive
                                dynamicmodel.lb(alpha) = 0;
                                dynamicmodel.ub(alpha) = 100;
                                dynamicmodel.obj(alpha) = -1*kappa2;
                                dynamicmodel.c(alpha) = -1*kappa2;
                                dynamicmodel.rxns{alpha} = 'alpha';
                                
                                % constraint 2
                                %     xi - ri <= eps2
                                %     ri >= 0
                                
                                % new row and column
                                rxnpos = size(dynamicmodel.A, 1) + 1;
                                alpha = size(dynamicmodel.A, 2) + 1;
                                dynamicmodel.A(rxnpos, rxn) = 1;
                                dynamicmodel.A(rxnpos, alpha) = -1;
                                dynamicmodel.rhs(rxnpos) = 0;
                                dynamicmodel.sense(rxnpos) = '<';
                                
                                % set ri to be positive
                                dynamicmodel.lb(alpha) = 0;
                                dynamicmodel.ub(alpha) = 100;
                                dynamicmodel.obj(alpha) = -1*kappa2;
                                dynamicmodel.rxns{alpha} = 'alpha';
                        end
                end
                dynamicmodel.vtype = repmat('C', size(dynamicmodel.A, 2), 1);
        end

        function [dynamicModelFlux] = computeDynamicModelSoln(dynamicmodel)
                values_to_query = [1:dynamicmodel.original_size];
                params.outputflag = 0;
                dynamicModelFlux = gurobi(dynamicmodel, params);
                dynamicModelFlux = dynamicModelFlux.x(values_to_query);
        end

gurobi_model = make_gurobi_model(model);
[time_points, model_positions, final_metabolomics] = ...
        preprocess_metabolomics(metabolomics_sheet, sheetname, timeseries);
slope = compute_fac(final_metabolomics, time_points);
normedSlope = normalize_fac(slope, normalizationFlag);
constrain_model = addMetabolomicConstraints(gurobi_model,...
        final_metabolomics, model_positions, normedSlope, kappa);
dynamicmodel = setDefaultConstraints(constrain_model, kappa2);
dynamicModelFlux = computeDynamicModelSoln(dynamicmodel);

%clear time_points model_positions final_metabolomics
end


