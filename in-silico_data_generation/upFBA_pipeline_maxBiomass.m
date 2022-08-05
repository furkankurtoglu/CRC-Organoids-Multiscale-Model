function [Model, stat, tmp_sol, v, r, p, q, Right_Mod] = upFBA_pipeline_maxBiomass(Model, initvalue, met_IDs, foldchange_means, foldchange_sds, SPECIES_BOUND)

    %% estimate metabolite rates of change
    changeSlopes = zeros(length(met_IDs), 1);
    changeIntervals = zeros(length(met_IDs), 1);
    SPECIES_HARD_CONSTRAINT = [];
    SPECIES_COUPLED_CONSTRAINT = [];
    for i = 1:length(met_IDs)
        tmp1 = initvalue(i) * (foldchange_means{i, 1} - 1) / 24.0;
        changeSlopes(i, 1) = tmp1;
        changeIntervals(i, 1) = 0;
        tmp_met_ID = split(met_IDs{i}, ";");
        if numel(tmp_met_ID) == 1
            SPECIES_HARD_CONSTRAINT = [SPECIES_HARD_CONSTRAINT; met_IDs(i) tmp1];
        else
            SPECIES_COUPLED_CONSTRAINT = [SPECIES_COUPLED_CONSTRAINT; met_IDs(i) tmp1];
        end
    end
    
    measurable_mets = [];
    for i = 1:length(met_IDs)
        measurable_mets = [measurable_mets; split(met_IDs{i}, "; ")];
    end
    
    %% change species bounds
    b_nonzero_indices_mat = zeros(89);  % 89 metabolites
    knockdown_biomass_mat = zeros(1); % 1 run (one set of initial conditions)
    
    Model = addSpeciesBounds_new(Model, SPECIES_BOUND);
    Model = addSpeciesHardConstraints_new(Model, SPECIES_HARD_CONSTRAINT);
    % optimize the model
%     for j=1:numel(model.rxns)
%         tmp_model = changeRxnBounds(model, model.rxns{j}, 0.0, 'u');
%         tmp_model = changeRxnBounds(tmp_model, tmp_model.rxns{j}, 0.0, 'l');
        tmp_model = Model;
        tmp_sol = optimizeCbModel(tmp_model, 'max');
        tmp_sol;

        save('tmp_sol.mat',"tmp_sol");

        if ~tmp_sol.stat
            relaxOption.internalRelax = 0;
            relaxOption.exchangeRelax = 0;
            relaxOption.steadyStateRelax = 1;
            relaxOption.epsilon = 10e-6;
            mets_to_exclude = false(size(tmp_model.mets));          
            mets_to_exclude(contains(tmp_model.mets, "_greater")) = true;
            mets_to_exclude(contains(tmp_model.mets, "_lower")) = true;
            relaxOption.excludedMetabolites = mets_to_exclude;
            solution = relaxedFBA(tmp_model,relaxOption);
            solution;
            [stat,v,r,p,q] = deal(solution.stat, solution.v, solution.r, solution.p, solution.q);
            if stat == 0
                v = zeros(size(tmp_model.rxns));
                r = zeros(size(tmp_model.mets));
                p = zeros(size(tmp_model.rxns));
                q = zeros(size(tmp_model.rxns));
                knockdown_biomass_mat(1) = NaN;
            else
                b_nonzero_indices = find(r ~= 0);
                b_nonzero_values = r(b_nonzero_indices);
                b_nonzero_indices_mat(b_nonzero_indices) = b_nonzero_indices_mat(b_nonzero_indices)+1;
                for m = 1:numel(b_nonzero_indices)
                    tmp_model.b(b_nonzero_indices(m)) = tmp_model.b(b_nonzero_indices(m)) - b_nonzero_values(m);
                end    
                solutionDel = optimizeCbModel(tmp_model, 'max');
                Right_Mod = tmp_model;
                stat = solutionDel.stat;
                if stat == 1
                   v = solutionDel.x;
                   if solutionDel.f > 0
                       knockdown_biomass_mat(1) = solutionDel.f;
                       %fprintf('m = %f\n', solutionDel.f)
                   else
                       knockdown_biomass_mat(1) = 0;
                   end
                else
                   v = zeros(size(tmp_model.rxns));
                   r = zeros(size(tmp_model.mets));
                   p = zeros(size(tmp_model.rxns));
                   q = zeros(size(tmp_model.rxns));
                   knockdown_biomass_mat(1) = NaN;
                end
            end
        else
            
            stat = tmp_sol.stat;
            %As there is no relaxation, no need to save relaxed FBA outputs
            v = zeros(size(tmp_model.rxns)); 
            r = zeros(size(tmp_model.mets));
            p = zeros(size(tmp_model.rxns));
            q = zeros(size(tmp_model.rxns));
            Right_Mod = tmp_model;
            if tmp_sol.f > 0
                knockdown_biomass_mat(1) = tmp_sol.f;
            else
                knockdown_biomass_mat(1) = 0;
            end
%             tmp_sol = ;
        end
%     end
    
end
