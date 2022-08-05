
function model = addSpeciesHardConstraints_new(model, SPECIES_HARD_CONSTRAINT)
    % obtain full stoichiometric matrix
    full_Model_S = full(model.S);

    % prepare data to be appended to Model
    species_S_rows = [];
    species_bs = [];
    species_mets = cell(size(SPECIES_HARD_CONSTRAINT, 1), 1);
    species_csenses = [];
    species_metNames = cell(size(SPECIES_HARD_CONSTRAINT, 1), 1);
    mets_to_remove = cell(size(SPECIES_HARD_CONSTRAINT, 1), 1);

    for i=1:size(SPECIES_HARD_CONSTRAINT, 1)
        species_metNames{i} = SPECIES_HARD_CONSTRAINT{i, 1};
        species_bs(i, 1) = SPECIES_HARD_CONSTRAINT{i, 2};
        species_csenses(i, 1) = 'E';
        species_mets{i} = SPECIES_HARD_CONSTRAINT{i, 1};
        [~, tmp_indices_bool] = ismember(model.mets, SPECIES_HARD_CONSTRAINT{i, 1});
        tmp_indices = find(tmp_indices_bool);
        species_S_rows(i, :) = full_Model_S(tmp_indices, :);
        mets_to_remove{i} = convertStringsToChars(SPECIES_HARD_CONSTRAINT{i, 1});
    end    

    % remove nonsteady state metabolites
    model = removeMetabolites(model, mets_to_remove, false);

    % append the newly generated data to Model
    model.S = sparse([full(model.S); species_S_rows]);
    model.b = [model.b; species_bs];
    model.mets = [model.mets; species_mets];
    model.csense = [model.csense; species_csenses];
    model.metNames = [model.metNames; species_metNames];
end
