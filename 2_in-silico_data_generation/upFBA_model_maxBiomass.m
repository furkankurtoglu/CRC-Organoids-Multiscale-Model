function [right] = upFBA_model_maxBiomass(init_value_data, fold_change_kras_data, fold_change_wt_data, folder_path, conditions, title_name, rate_vals)

%     HK_kd = 'N';


    num_of_met_run = 5;
    save_data = "Y";
    Simulate_WT = "Y";
    Simulate_KRAS = "N";
    LHS = "N";


    %% write reactions
    reactionFormulas = {'M_glc_D_e -> M_glc_D_c',...  % GLUT
    'M_glc_D_c + M_atp_i -> M_g6p_c + M_adp_i + M_h_c',...  % HK
    'M_g6p_c <=> M_f6p_c',...  % HPI
    'M_f6p_c + M_atp_i -> M_fdp_c + M_adp_i + M_h_c',...  % PFK
    'M_fdp_c + M_h2o_c -> M_f6p_c + M_pi_c',...  % FBP (gluconeogenesis)
    'M_fdp_c <=> M_dhap_c + M_g3p_c',...  % ALDO
    'M_dhap_c <=> M_g3p_c',...  % TPI
    'M_g3p_c + M_nad_i + M_pi_c <=> M_13dpg_c + M_h_c + M_nadh_c',...  % GAPDH
    'M_13dpg_c + M_h2o_c -> M_3pg_c + M_h_c + M_pi_c',...  % 13DPG-phosphatase 
    'M_13dpg_c + M_adp_i <=> M_3pg_c + M_atp_i',...  % PGK
    'M_13dpg_c <=> M_23dpg_c + M_h_c',...  % DPGM (mutase, Rapoport?Luebering-pathway)
    'M_23dpg_c + M_h2o_c -> M_3pg_c + M_pi_c',...  % DPGase
    'M_3pg_c <=> M_2pg_c',...  % PGAM
    'M_2pg_c <=> M_h2o_c + M_pep_c',...  % ENO
    'M_pep_c + M_adp_i + M_h_c -> M_pyr_c + M_atp_i',...  % PYK
    'M_pyr_m + M_coa_m + M_nad_i -> M_accoa_m + M_co2_m + M_nadh_m',...  % PDH (pyruvate enters mitochondria)
    'M_g6p_c + M_nadp_c <=> M_6pgl_c + M_h_c + M_nadph_c',...  % G6PDH
    'M_6pgl_c + M_h2o_c -> M_6pgc_c + M_h_c',...  % PGL (lactonase)
    'M_6pgc_c + M_nadp_c -> M_r5p_ru5p_c + M_nadph_c + M_co2_c',...  % 6PGDH
    'M_r5p_ru5p_c <=> M_xu5p_D_c',...  % RPE
    'M_r1p_c <=> M_r5p_ru5p_c',...  % RPI
    'M_g3p_c + M_s7p_c <=> M_f6p_c + M_e4p_c',...  % TA
    'M_r5p_ru5p_c + M_xu5p_D_c <=> M_g3p_c + M_s7p_c',...  % TK1
    'M_xu5p_D_c + M_e4p_c <=> M_f6p_c + M_g3p_c',...  % TK2
    'M_accoa_m + M_oaa_i + M_h2o_m -> M_cit_icit_i + M_h_m + M_coa_m',...  % CS (citrate synthase)
    'M_atp_i + M_coa_c + M_cit_icit_i -> M_adp_i + M_pi_c + M_accoa_c + M_oaa_i',...  % CLY (citrate lyase)
    'M_nadp_c + M_cit_icit_i -> M_nadph_c + M_akg_i + M_co2_c',...  % IDH1 
    'M_nad_i + M_cit_icit_i -> M_akg_i + M_co2_m + M_nadh_m',...  % IDH3
    'M_nadp_m + M_cit_icit_i <=> M_nadph_m + M_akg_i + M_co2_m',...  % IDH2 
    'M_akg_i + M_coa_m + M_nad_i -> M_succoa_m + M_co2_m + M_nadh_m',...  % AKGD
    'M_adp_i + M_pi_m + M_succoa_m <=> M_coa_m + M_atp_i + M_succ_m',...  % SCOAS
    'M_fad_m + M_succ_m <=> M_fadh2_m + M_fum_i',...  % SDHA
    'M_q10_m + M_succ_m -> M_q10h2_m + M_fum_i',...  % SDHC+SDHD
    'M_fum_i + M_h2o_c <=> M_mal_L_i',...  % FUM2 
    'M_fum_i + M_h2o_m <=> M_mal_L_i',...  % FUM1 
    'M_mal_L_i + M_nad_i <=> M_oaa_i + M_h_c + M_nadh_c',...  % MDH2 
    'M_mal_L_i + M_nad_i <=> M_oaa_i + M_h_m + M_nadh_m',...  % MDH1 
    'M_mal_L_i + M_nadp_m <=> M_pyr_m + M_co2_m + M_h_m + M_nadph_m',...  % MMALIC
    'M_fadh2_m + M_q10_m -> M_fad_m + M_q10h2_m',...  % succinate-Q reductase
    '2 M_h_m + 0.5 M_o2_m + M_q10h2_m -> M_h2o_m + M_q10_m + 2 M_h_om',...  % cytochrome oxidase
    '5 M_h_m + M_nadh_m + M_q10_m -> 4 M_h_om + M_nad_i + M_q10h2_m',...  % NADH-Q reductase
    'M_pi_m + M_adp_i + 4 M_h_om -> M_h2o_m + 3 M_h_m + M_atp_i',...  % ATP synthase
    'M_pyr_c + M_h_c + M_nadh_c <=> M_lac_L_c + M_nad_i',...  % LDH
    'M_h_c + M_lac_L_c <=> M_h_e + M_lac_L_e',...  % MCT
    'M_h2o_m + M_gln_L_i -> M_glu_L_m + M_nh4_m',...  % GLNS (glutaminase)
    'M_h2o_m + M_nad_i + M_glu_L_m <=> M_h_m + M_akg_i + M_nadh_m + M_nh4_m',...  % GDH2
    'M_h2o_m + M_nadp_m + M_glu_L_m <=> M_h_m + M_akg_i + M_nadph_m + M_nh4_m',...  % GDH1+GDH3
    'M_h2o_c <=> M_h2o_m',...  % h2o m <=> c
    'M_o2_c <=> M_o2_m',...  % o2 m <=> c
    'M_co2_c <=> M_co2_m',...  % co2 m <=> c
    'M_pi_c <=> M_pi_m',...  % pi m <=> c
    'M_h_c <=> M_h_m',...  % h m <=> c
    'M_pyr_c <=> M_pyr_m',...  % MPC
    'M_nadh_c <=> M_nadh_m',...  % nadh m <=> c
    'M_nh4_c <=> M_nh4_m',...  % nh4 m <=> c
    'M_h2o_e <=> M_h2o_c',...  % H2O transport
    'M_o2_e <=> M_o2_c',...  % O2 transport
    'M_co2_e <=> M_co2_c',...  % CO2 transport
    'M_pi_e <=> M_pi_c',...  % phosphate transport
    'M_h_e <=> M_h_c',...  % H+ transport
    'M_gln_L_e <=> M_gln_L_i',...  % ASCT2 (Gln transport)
    'M_nh4_e <=> M_nh4_c',...  % NH4 transport
    'M_h2o_b <=> M_h2o_e',...  % H2O exchange
    'M_o2_b <=> M_o2_e',...  % O2 exchange
    'M_co2_b <=> M_co2_e',...  % CO2 exchange
    'M_pi_b <=> M_pi_e',...  % phosphate exchange
    'M_h_b <=> M_h_e',...  % H+ exchange
    'M_lac_L_b <=> M_lac_L_e',...  % lactate exchange
    'M_glc_D_b <=> M_glc_D_e',...  % glucose exchange
    'M_gln_L_b <=> M_gln_L_e',...  % glutamine exchange
    'M_nh4_b <=> M_nh4_e',...  % nh4 exchange
    'M_atp_i + M_h2o_c -> M_adp_i + M_h_c + M_pi_c',...  % ATPase
    '0.87841 M_3pg_c + 0.536023 M_pyr_c + 3.425741 M_oaa_i + 1.62 M_akg_i + 3.17174 M_accoa_m + 0.16035 M_r5p_ru5p_c + 100 M_atp_i + 100 M_h2o_c + 18 M_nadph_c  + 3.5 M_nad_i -> biomass + 3.17174 M_coa_m + 100 M_adp_i + 100 M_h_c + 100 M_pi_c + 18 M_nadp_c + 3.5 M_nadh_c',...  % biomass
    'M_glu_L_m + M_atp_i + M_nh4_m -> M_gln_L_i + M_h_m + M_adp_i + M_pi_m'...  % glutamine synthetase (GS)
    }; 

    reactionNames = {'GLUT', 'HK', 'HPI', 'PFK', 'FBP', 'ALDO', ...
        'TPI', 'GAPDH', '13DPG-phophatase', 'PGK', 'DPGM', 'DPGase', ...
        'PGAM', 'ENO', 'PYK', 'PDH', 'G6PDH', 'PGL', ...
        '6PGDH', 'RPE', 'RPM', 'TA', 'TK1', 'TK2', ...
        'CS', 'CLY', 'IDH1', 'IDH3', 'IDH2', 'AKGD', ...
        'SCOAS', 'SDHA', 'SDHC(D)', 'FUM2', 'FUM1', 'MDH2', ...
        'MDH1', 'MMALIC', 'succinate-Q reductase', 'cytochrome oxidase', 'NADH-Q reductase', 'ATP synthase', ...
        'LDH', 'MCT', 'GLNS', 'GDH2', 'GDH1(3)', 'H2O (m-c)', ...
        'O2 (m-c)', 'CO2 (m-c)', 'Pi (m-c)', 'H (m-c)', 'MPC', 'NADH (m-c)', ...
        'NH4 (m-c)', 'H2O (e-c)', 'O2 (e-c)', 'CO2 (e-c)', 'Pi (e-c)', 'H (e-c)', ...
        'ASCT2', 'NH4 (e-c)', 'H2O (b-e)', 'O2 (b-e)', 'CO2 (b-e)', 'Pi (b-e)', ...
        'H (b-e)', 'Lac (b-e)', 'Glc (b-e)', 'Gln (b-e)', 'NH4 (b-e)', 'ATPase', ...
        'biomass_obj', 'GS'}; 

    % specify bounds for flux
    lowerbounds = -500 * ones(size(reactionNames));
    lowerbounds([1, 2, 4, 5, 9, 12, ...
        15, 16, 18, 19, 25, 26, ...
        27, 28, 30, 33, ...
        39, 40, 41, 42, 45, ...
        72, 73, 74]) = 0; 
    upperbounds = 500 * ones(size(reactionNames));

%     lowerbounds(2) = 0;
%     upperbounds(2) = 0;

%     lowerbounds(1) = rate_vals(1); % GLUT
%     lowerbounds(44) = rate_vals(2); % MCT
%     lowerbounds(61) = rate_vals(3); % ASNT2
%     upperbounds(1) = rate_vals(1); % GLUT
%     upperbounds(44) = rate_vals(2); % MCT
%     upperbounds(61) = rate_vals(3); % ASNT2

    % construct the model (preliminary)
    Model = createModel(reactionNames, reactionNames, reactionFormulas, ...
        'lowerBoundList', lowerbounds, 'upperBoundList', upperbounds);
    Model = addCOBRAConstraints(Model, {'HPI', 'G6PDH'}, 0, 'c', ...
        [1, -9], 'dsense', ['E']);
    Model = changeObjective(Model, 'biomass_obj');



    %% load experimental data
    met_IDs_wt = fold_change_wt_data{:, 1};
    met_IDs_kras = fold_change_kras_data{:, 1};
    foldchange_means_wt = fold_change_wt_data(:, 2);
    foldchange_means_kras = fold_change_kras_data(:, 2);
    foldchange_sds_wt = fold_change_wt_data(:, 3);
    foldchange_sds_kras = fold_change_kras_data(:, 3);
    
    %% clear all unsteady-state assumptions
    %% perform upFBA ------ Glycolysis, PPP, and TCA cycle
    


    num_of_FBA_runs = num_of_met_run*num_of_met_run;
    WT_Model = struct;
    WT_Model.model_lst = cell(1, num_of_FBA_runs);
    WT_Model.stat_lst = zeros(1, num_of_FBA_runs);
    WT_Model.sol_lst = cell(1, num_of_FBA_runs);
    WT_Model.v_lst = [];
    WT_Model.r_lst = [];
    WT_Model.p_lst = [];
    WT_Model.q_lst = [];
    WT_Model.initvalue_lst = [];
    
    %% generate random initial values
    %init_value_mat = init_value_data(:, 1) + (init_value_data(:, 2) - init_value_data(:, 1)) ...
       %     .* rand(size(init_value_data(:, 1), 1), num_runs);



    
    lac_bounds = [6.4, 19.2];
    gln_bounds = [0.56, 1.08];
    glu_bounds = [0.88, 1.56];

    init_value_mat = zeros(7,num_of_FBA_runs);
    init_value_mat(1,:) = 3.87;
    init_value_mat(2,:) = 0.62;
    init_value_mat(3,:) = 0.54;
    init_value_mat(4,:) = 0.54;
    init_value_mat(5,:) = sum(lac_bounds)/2;
    init_value_mat(6,:) = sum(gln_bounds)/2;
    init_value_mat(7,:) = sum(glu_bounds)/2;

    if (LHS == 'N')

%         glc_bds = linspace(-0.223, 0,num_of_met_run);
%         gln_bds = linspace(-0.003, 0,num_of_met_run);
%     
%         g6p_i_bds = linspace(1.68,6.0,num_of_met_run);
%         fbp_i_bds = linspace(0.16,1.08,num_of_met_run);
%         g3p_i_bds = linspace(0.24,0.84,num_of_met_run);
%         pep_i_bds = linspace(0.24,0.84,num_of_met_run);
%         lac_i_bds = linspace(6.4,19.2,num_of_met_run);
%         gln_i_bds = linspace(0.56,1.08,num_of_met_run); 
%         glu_i_bds = linspace(0.88,1.56,num_of_met_run);


        glc_bds = linspace(-0.223, 0,num_of_met_run);
        gln_bds = linspace(-0.003, 0,num_of_met_run);
    
        g6p_i_bds = linspace(1.68,6.0,num_of_met_run);
        fbp_i_bds = linspace(0.16,1.08,num_of_met_run);
        g3p_i_bds = linspace(0.24,0.84,num_of_met_run);
        pep_i_bds = linspace(0.24,0.84,num_of_met_run);
        lac_i_bds = linspace(6.4,19.2,num_of_met_run);
        gln_i_bds = linspace(0.56,1.08,num_of_met_run); 
        glu_i_bds = linspace(0.88,1.56,num_of_met_run);


        
        %%delete('data.csv');
        save('init_vals.mat',"init_value_mat");

        if save_data == "Y"
            names = {'Glucose_exchange_input','Glutamine_exchange_input','G6P_input','FBP_input','G3P_input','PEP_input','lac_exchange_input','gln_input','glu_input','solution_output'};
%             names = {'Glucose_exchange_input','Glutamine_exchange_input','solution_output'};
            writecell(names,'../3_DNN_model_generation/WT_in_silico_data_GLC_GLN_and_Seven_Metabolites_5_even.csv',Delimiter=',',WriteMode='append')
        end
        i=1;
        total_iteration = 0;
        if (Simulate_WT == "Y")
            for k = 1:num_of_met_run %1
                for l = 1:num_of_met_run %2
                    for m = 1:num_of_met_run %3 
                       for n = 1:num_of_met_run %4
                           for o = 1:num_of_met_run %5
                               for p = 1:num_of_met_run %6
                                   for q = 1:num_of_met_run %7
                                        for r = 1:num_of_met_run %7
                                            for s = 1:num_of_met_run %8
                                               % for t = 1:num_of_met_run %10 (HK)
                                    % specify bounds for the metabolite concentrations
                                    SPECIES_BOUND_WT={
                                    "M_h2o_b[c]"	-500	500	;	% 1 M_h2o_b
                                    "M_o2_b[c]"	-500	500	;	% 2 M_o2_b
                                    "M_co2_b[c]"	-500	500	;	% 3 M_co2_b
                                    "M_pi_b[c]"	-500	500	;	% 4 M_pi_b
                                    "M_h_b[c]"	-500	500	;	% 5 M_h_b
                                    "M_lac_L_b[c]"	0.234	0.234;	% 6 M_lac_L_b
                                    "M_glc_D_b[c]"	glc_bds(k)	glc_bds(k)	;	% 7 M_glc_D_b
                                    "M_gln_L_b[c]"	gln_bds(l)	gln_bds(l);	% 8 M_gln_L_b
                                    "M_nh4_b[c]"	-500	500	;	% 9 M_nh4_b
                                    %"biomass[c]" rate_vals(7) rate_vals(7) ; % growth rate as constraint; baseline
                                    "biomass[c]" 0   500 ; % 10 biomass
                                };

                                initvalue = init_value_mat(:, 1);
                                initvalue(5) = g6p_i_bds(m);
                                initvalue(5) = fbp_i_bds(n);
                                initvalue(5) = g3p_i_bds(o);
                                initvalue(5) = pep_i_bds(p);
                                initvalue(5) = lac_i_bds(q);
                                initvalue(6) = gln_i_bds(r);
                                initvalue(7) = glu_i_bds(s); 

                                [new_Model, stat, sol, v1, r1, p1, q1, right,iteration] = upFBA_pipeline_maxBiomass(Model, initvalue, met_IDs_wt, foldchange_means_wt, foldchange_sds_wt, SPECIES_BOUND_WT);
                                %WT_Model.model_lst{i} = new_Model;
                                %WT_Model.stat_lst(i) = stat;
                                %WT_Model.sol_lst{i} = sol;
                                %size(sol.x)
                                if (size(sol.x)>0)
                                    solution = sol.x(73);
                                else
                                    solution = 0.0;
                                end
                                %WT_Model.v_lst = [WT_Model.v_lst v];
                                %WT_Model.r_lst = [WT_Model.r_lst r];
                                %WT_Model.p_lst = [WT_Model.p_lst p];
                                %WT_Model.q_lst= [WT_Model.q_lst q];
                                %WT_Model.initvalue_lst{i} = [WT_Model.initvalue_lst initvalue];
                                i=i+1;
                                %fprintf('m = %i\t %i\t %i\t %i\n', k,m,lac_sol,solution)
                                %sol_vector = [-1*glc_bds(k),-1*gln_bds(l),lac_i_bds(m),gln_i_bds(n),glu_i_bds(o),solution,intra_glu_change,intra_gln_change,intra_lac_change];
                                %sol_vector = [-1*glc_bds(k),-1*gln_bds(l),solution];
                                int_met_vector = [-1*glc_bds(k),-1*gln_bds(l),g6p_i_bds(m),fbp_i_bds(n),g3p_i_bds(o),pep_i_bds(p),lac_i_bds(q),gln_i_bds(r),glu_i_bds(s),solution];
                                %int_met_vector = [-1*glc_bds(k),-1*gln_bds(l),solution];
                                if save_data == "Y"
                                    %%writematrix(sol_vector,'./data.csv',Delimiter=',',WriteMode='append');
                                    %writematrix(sol_vector,'../DNN_model_generation/WT_in_silico_data_even_NO_INTRACELLULAR.csv',Delimiter=',',WriteMode='append');
                                    writematrix(int_met_vector,'../3_DNN_model_generation/WT_in_silico_data_GLC_GLN_and_Seven_Metabolites_5_even.csv',Delimiter=',',WriteMode='append');
                                    %writematrix(int_met_vector,'../3_DNN_model_generation/HK_Gene_Knockout_trials.csv',Delimiter=',',WriteMode='append')
                                end
                                total_iteration = total_iteration + iteration;
%                                                 end
                                            end
                                        end
                                   end
                               end
                            end
                        end
                    end
                end
%                 fprintf('m = %i\t %i\t %i\t %i\n', k,l,lac_sol,solution)
                %total_iteration
            end
        end
        total_iteration

    end
    i=1;

    KRAS_num_of_met_run = num_of_met_run;
    if (Simulate_KRAS == "Y")
        for k = 1:KRAS_num_of_met_run
            for l = 1:KRAS_num_of_met_run
                for m = 1:KRAS_num_of_met_run
                    for n = 1:KRAS_num_of_met_run
                        for o = 1:KRAS_num_of_met_run
                            SPECIES_BOUND_KRAS={
                            "M_h2o_b[c]"	-500	500	;	% 1 M_h2o_b
                            "M_o2_b[c]"	-500	500	;	% 2 M_o2_b
                            "M_co2_b[c]"	-500	500	;	% 3 M_co2_b
                            "M_pi_b[c]"	-500	500	;	% 4 M_pi_b
                            "M_h_b[c]"	-500	500	;	% 5 M_h_b
                            "M_lac_L_b[c]"	0.234	0.234;	% 6 M_lac_L_b
                            "M_glc_D_b[c]"	KRAS_glc_bds(k)	KRAS_glc_bds(k)	;	% 7 M_glc_D_b
                            "M_gln_L_b[c]"	KRAS_gln_bds(l)	KRAS_gln_bds(l);	% 8 M_gln_L_b
                            "M_nh4_b[c]"	-500	500	;	% 9 M_nh4_b
                            %"biomass[c]" rate_vals(7) rate_vals(7) ; % growth rate as constraint; baseline
                            "biomass[c]" 0   500;  % 10 biomass
                             };
                            initvalue = init_value_mat(:, 1);
                            initvalue(5)= KRAS_lac_i_bds(m);
                            initvalue(6)= KRAS_gln_i_bds(n);
                            initvalue(7)= KRAS_glc_i_bds(o);
                            [new_Model, stat, sol, v, r, p, q] = upFBA_pipeline_maxBiomass(Model_2, initvalue, met_IDs_kras, foldchange_means_kras, foldchange_sds_kras, SPECIES_BOUND_KRAS);
%                             KRAS_Model.model_lst{i} = new_Model;
%                             KRAS_Model.stat_lst(i) = stat;
%                             KRAS_Model.sol_lst{i} = sol;
%                             KRAS_Model.v_lst = [KRAS_Model.v_lst v];
%                             KRAS_Model.r_lst = [KRAS_Model.r_lst r];
%                             KRAS_Model.p_lst = [KRAS_Model.p_lst p];
%                             KRAS_Model.q_lst = [KRAS_Model.q_lst q];  
%                             KRAS_Model.initvalue_lst = [KRAS_Model.initvalue_lst initvalue];
                            if (size(sol.x)>0)
                                lac_sol = sol.x(68);
                                solution = sol.x(73);
                                intra_lac_change = sol.x(44)-sol.x(43);
                                intra_g6p_change = sol.x(1)-sol.x(2);
                                intra_gln_change = sol.x(61)+sol.x(74)-sol.x(45);
                            else
                                solution = 0.0;
                                intra_lac_change = 0;
                                intra_g6p_change = 0;
                                intra_gln_change = 0;
                            end
                            i=i+1;
                           % fprintf('m = %i\t %i\t %i\t %i\n', k,m,lac_sol,solution)
                            sol_vector = [-1*KRAS_glc_bds(k),-1*KRAS_gln_bds(l),KRAS_lac_i_bds(m),KRAS_gln_i_bds(n),KRAS_glc_i_bds(o),solution,intra_g6p_change,intra_gln_change,intra_lac_change];
                            if save_data == "Y"
                                writematrix(sol_vector,'../DNN_model_generation/KRAS_in_silico_data.csv',Delimiter=',',WriteMode='append');
                            end
                        end
                    end
                end
            end
            fprintf('m = %i\t %i\t %i\t %i\n', k,l,lac_sol,solution)
        end
        fclose all;
end
