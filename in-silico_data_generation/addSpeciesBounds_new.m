
function model = addSpeciesBounds_new(model, SPECIES_BOUND)
    % obtain full stoichiometric matrix
    full_Model_S = full(model.S);

    % prepare data to be appended to Model
    species_lb_S_rows = [];
    species_ub_S_rows = [];
    species_lb_bs = [];
    species_ub_bs = [];
    species_lb_mets = cell(size(SPECIES_BOUND, 1), 1);
    species_ub_mets = cell(size(SPECIES_BOUND, 1), 1);
    species_lb_csenses = [];
    species_ub_csenses = [];
    species_lb_metNames = cell(size(SPECIES_BOUND, 1), 1);
    species_ub_metNames = cell(size(SPECIES_BOUND, 1), 1);
    mets_to_remove = cell(size(SPECIES_BOUND, 1), 1);

    for i=1:size(SPECIES_BOUND, 1)
        species_lb_metNames{i} = convertStringsToChars(strcat(SPECIES_BOUND{i, 1}, '_greater'));
        species_ub_metNames{i} = convertStringsToChars(strcat(SPECIES_BOUND{i, 1}, '_lower'));
        species_lb_bs(i, 1) = SPECIES_BOUND{i, 2};
        species_ub_bs(i, 1) = SPECIES_BOUND{i, 3};
        species_lb_csenses(i, 1) = 'G';
        species_ub_csenses(i, 1) = 'L';
        species_lb_mets{i} = convertStringsToChars(strcat(SPECIES_BOUND{i, 1}, '_greater'));
        species_ub_mets{i} = convertStringsToChars(strcat(SPECIES_BOUND{i, 1}, '_lower'));
        [~, tmp_indices_bool] = ismember(model.mets, SPECIES_BOUND{i, 1});
        tmp_indices = find(tmp_indices_bool);
        species_lb_S_rows(i, :) = full_Model_S(tmp_indices, :);
        species_ub_S_rows(i, :) = full_Model_S(tmp_indices, :);
        mets_to_remove{i} = convertStringsToChars(SPECIES_BOUND{i, 1});
    end    

    % remove nonsteady state metabolites    
    model = removeMetabolites(model, mets_to_remove, false);

    % append the newly generated data to Model
    model.S = sparse([full(model.S); species_lb_S_rows; species_ub_S_rows]);
    model.b = [model.b; species_lb_bs; species_ub_bs];
    model.mets = [model.mets; species_lb_mets; species_ub_mets];
    model.csense = [model.csense; species_lb_csenses; species_ub_csenses];
    model.metNames = [model.metNames; species_lb_metNames; species_ub_metNames];
end
