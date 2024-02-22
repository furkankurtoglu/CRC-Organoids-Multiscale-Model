clear all
close all
clc

changeCobraSolver('glpk');

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

SPECIES_BOUND_WT={
                    "M_h2o_b[c]"	-500	500	;	% 1 M_h2o_b
                    "M_o2_b[c]"	-500	500	;	% 2 M_o2_b
                    "M_co2_b[c]"	-500	500	;	% 3 M_co2_b
                    "M_pi_b[c]"	-500	500	;	% 4 M_pi_b
                    "M_h_b[c]"	-500	500	;	% 5 M_h_b
                    "M_lac_L_b[c]"	0.283	0.283;	% 6 M_lac_L_b
                    "M_glc_D_b[c]"	0.223	0.223	;	% 7 M_glc_D_b
                    "M_gln_L_b[c]"	0.003	0.003;	% 8 M_gln_L_b
                    "M_nh4_b[c]"	-500	500	;	% 9 M_nh4_b
                    "biomass[c]" 0   500 ; % 10 biomass
                };


Model = addSpeciesBounds_new(Model, SPECIES_BOUND_WT);
FBA_Sol = optimizeCbModel(Model);