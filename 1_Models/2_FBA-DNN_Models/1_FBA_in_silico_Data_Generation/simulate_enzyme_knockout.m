

initCobraToolbox(false)

folder_path_prefix = '.';  % specify the location of "data.mat" (flux distributions of the original phenotype)
folder_path_lst = {[folder_path_prefix 'WT_vs_KRAS_CRC_media/'], [folder_path_prefix 'WT_vs_KRAS_CAF_media/'], ...
    [folder_path_prefix 'WT_CRC_media_vs_CAF_media/'], [folder_path_prefix 'KRAS_CRC_media_vs_CAF_media/']};  % these four folders should have been created. 
conditions_lst = {{'WT', 'KRAS'}, {'WT', 'KRAS'}, ...
    {'CRC media', 'CAF media'}, {'CRC media', 'CAF media'}};
title_name_lst = {'WT vs KRAS in CRC media', 'WT vs KRAS in CAF media', ...
    'WT in CRC media vs CAF media', 'KRAS in CRC media vs CAF media'};

b_nonzero_indices_mat = zeros(4, 100, 74, 89);  % 4 conditions, 100 sets of mass balance constraints, 74 reactions, 89 metabolites
knockdown_biomass_mat = zeros(4, 100, 74);
for i=1:4    % i=1, WT(CRC-only media); i=2, KRAS(CRC-CAF media); i=3, WT(CRC-CAF media); i=4, KRAS(CRC-only media)
    data_path = [folder_path_lst{i} 'data.mat'];
    load(data_path);
    for k=1:100
        disp(k)
        if (i == 2) || (i == 3)
            model = KRAS_Model.model_lst{1, k};
            model.b(82) = KRAS_Model.v_lst(73, k);
        else
            model = WT_Model.model_lst{1, k};
            model.b(82) = WT_Model.v_lst(73, k);
        end
        model.b(72) = 0;
        for j=1:numel(model.rxns)
            tmp_model = changeRxnBounds(model, model.rxns{j}, 0.0, 'u');
            tmp_model = changeRxnBounds(tmp_model, tmp_model.rxns{j}, 0.0, 'l');
            tmp_sol = optimizeCbModel(tmp_model, 'max');
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
                [stat,v,r,p,q] = deal(solution.stat, solution.v, solution.r, solution.p, solution.q);
                if stat == 0
                    v = zeros(size(tmp_model.rxns));
                    r = zeros(size(tmp_model.mets));
                    p = zeros(size(tmp_model.rxns));
                    q = zeros(size(tmp_model.rxns));
                    knockdown_biomass_mat(i, k, j) = NaN;
                else
                    b_nonzero_indices = find(r ~= 0);
                    b_nonzero_values = r(b_nonzero_indices);
                    b_nonzero_indices_mat(i, k, j, b_nonzero_indices) = b_nonzero_indices_mat(i, k, j, b_nonzero_indices)+1;
                    for m = 1:numel(b_nonzero_indices)
                        tmp_model.b(b_nonzero_indices(m)) = tmp_model.b(b_nonzero_indices(m)) - b_nonzero_values(m);
                    end    
                    solutionDel = optimizeCbModel(tmp_model, 'max');
                    stat = solutionDel.stat;
                    if stat == 1
                       v = solutionDel.x;                      
                       if solutionDel.f > 0
                           knockdown_biomass_mat(i, k, j) = solutionDel.f;
                       else
                           knockdown_biomass_mat(i, k, j) = 0;
                       end
                    else
                       v = zeros(size(tmp_model.rxns));
                       r = zeros(size(tmp_model.mets));
                       p = zeros(size(tmp_model.rxns));
                       q = zeros(size(tmp_model.rxns));
                       knockdown_biomass_mat(i, k, j) = NaN;
                    end
                end
            else
                if tmp_sol.f > 0
                    knockdown_biomass_mat(i, k, j) = tmp_sol.f;
                else
                    knockdown_biomass_mat(i, k, j) = 0;
                end
            end
        end
    end
end

save("enzyme_knockout_result.mat", "knockdown_biomass_mat", "b_nonzero_indices_mat");

