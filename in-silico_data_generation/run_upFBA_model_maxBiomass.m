
init_value_data_full = readtable('measured_initial_concentrations_boundaries.csv');
init_value_data = table2array(init_value_data_full(:, [2 3]));

fold_change_wt_data_lst = {'fold_change_wt_crc.csv', ...
    'fold_change_wt_caf.csv', ...
    'fold_change_wt_crc.csv', ... 
    'fold_change_kras_crc.csv'};
fold_change_kras_data_lst = {'fold_change_kras_crc.csv', ...
    'fold_change_kras_caf.csv', ...
    'fold_change_wt_caf.csv', ...
    'fold_change_kras_caf.csv'};
folder_path_prefix = '';   % set where you want to place the output files
folder_path_lst = {[folder_path_prefix 'max_WT_vs_KRAS_CRC_media/'], [folder_path_prefix 'max_WT_vs_KRAS_CAF_media/'], ...
    [folder_path_prefix 'max_WT_CRC_media_vs_CAF_media/'], [folder_path_prefix 'max_KRAS_CRC_media_vs_CAF_media/']};  % create these four directories before running the codes below
conditions_lst = {{'WT', 'KRAS'}, {'WT', 'KRAS'}, ...
    {'CRC media', 'CAF media'}, {'CRC media', 'CAF media'}};
title_name_lst = {'WT vs KRAS in CRC media', 'WT vs KRAS in CAF media', ...
    'WT in CRC media vs CAF media', 'KRAS in CRC media vs CAF media'};

% glucose, lactate, and glutamine uptake/secretion rates; biomass growth rates
rate_vals_mat = [0.223 0.283 0.003 0.210 0.234 0.003 0.035 0.034;  ... % WT(CRC) vs KRAS(CRC)
    0.573 0.784 0.052 0.927 1.117 0.105 0.033 0.033; ... % WT(CAF) vs KRAS(CAF)
    0.223 0.283 0.003 0.573 0.784 0.052 0.035 0.033; ... % WT(CRC) vs WT(CAF)
    0.210 0.234 0.003 0.927 1.117 0.105 0.034 0.033];  % KRAS(CRC) vs KRAS(CAF)

    i=1;
    rate_vals = rate_vals_mat(i, :);
    fold_change_wt_data = readtable(fold_change_wt_data_lst{i});
    fold_change_kras_data = readtable(fold_change_kras_data_lst{i});
    folder_path = folder_path_lst{i};
    conditions = conditions_lst{i};
    title_name = title_name_lst{i};
    WT = upFBA_model_maxBiomass(init_value_data, fold_change_kras_data, ...
        fold_change_wt_data, folder_path, conditions, title_name, rate_vals);
