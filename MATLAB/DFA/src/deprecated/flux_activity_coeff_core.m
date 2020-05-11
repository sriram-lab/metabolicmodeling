function [dynamicmodel,geneko_growthrate_obj, geneko_growthrate,core_values,manual_matchpos,headers] = flux_activity_coeff_core(model, timecourse_metabolomics_datafile,celllines)

% for CORE values - extracellular metabolomics over time.

model_rpmi = model;

% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % change objective if necessary to atp.
% model_rpmi = addExchangeRxn(model_rpmi,{'atp[c]'},0,100); %
% model_rpmi.c(find(model_rpmi.c)) = 0;
% model_rpmi.c(end) = 1;
% % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% load metabolomics data
[num txt] = xlsread(timecourse_metabolomics_datafile);
manual_matchpos = num(:,1); % coresponding position in model
core_values = num(:,2:end);
headers = txt(1,3:end)';

ix0 = (manual_matchpos(:,1) ~= 0);
tol1 = 0.01;

if (~exist('celllines','var')) || (isempty(celllines))
    celllines = headers;
end

kappa = 1;



clear model
model1 = model_rpmi;
model = model1;
model.A = model1.S;
model.obj = model1.c;
model.rhs = model1.b;
model.sense =repmat( '=',[size(model1.S,1),1]);
model.lb = model1.lb;
model.ub = model1.ub;
model.vtype = repmat('C',size(model1.S,2),1);
model.modelsense = 'max';

unqgenes = unique(model.genes);
for kk = 1:length(unqgenes),
    model2 = deleteModelGenes(model,unqgenes(kk));
    
    for j = 1:size(core_values,2)
        if ismember(headers(j), celllines)
            
            m2 = model2;%_rpmi;
            for i = 1:length(core_values(:,j))
                if ix0(i)
                    u3pos = manual_matchpos(i,1);
                    tol11(i) = core_values(i,j)/max(abs(core_values(:))); % FLUX ACTIVITY COEFFICIENT. normalized.
                    rowpos = u3pos;
                    colpos = size(m2.A,2) + 1;
                    colpos1 = colpos + 1;
                    
                    epsilon2 = tol11(i);
                    m2.rhs(rowpos) = epsilon2;
                    
                    m2.A(rowpos,colpos) = 1;
                    m2.A(rowpos,colpos1) = -1;
                    m2.vtype(colpos1) = 'C';
                    m2.vtype(colpos) = 'C';
                    
                    m2.lb(colpos1) = 0;
                    m2.ub(colpos1) = 1000;
                    m2.lb(colpos) = 0;
                    m2.ub(colpos) = 1000;
                    
                    m2.obj(colpos) = -1*kappa; % minimized
                    m2.obj(colpos1) = -1*kappa; % minimized
                    
                end
            end
            
            m2.vtype = repmat('C',size(m2.A,2),1);
            params.outputflag = 0;
            
            dynamicmodel = m2;
            
            solg1 = gurobi(m2,params);
            geneko_growthrate(kk,j) = solg1.x(find(model_rpmi.c));
            geneko_growthrate_obj(kk,j) = solg1.objval;
            disp(kk)

        end
        
        
    end
end

end