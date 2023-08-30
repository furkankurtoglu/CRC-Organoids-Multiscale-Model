%Net Export rate calculations
clear all
close
clc


glc_con = 17.5; %mM
gln_con = 5.5; %mM


time = 24; %hr
integrated_cell_count = 440000; %cells ???

volume_of_well = 1.5e-3; %l



glucose_uptake_rate = 0.223; %mM/hr = mmol/l/hr
glutamine_uptake_rate = 0.003;



e_0_glu_wt_crc = glucose_uptake_rate * volume_of_well / glc_con / integrated_cell_count;
e_0_gln_wt_crc = glutamine_uptake_rate * volume_of_well / gln_con /integrated_cell_count;


glc_cons = 0:0.1:glucose_uptake_rate;

uptake_


% uptaken_glucose_conc_by_rate = glucose_uptake_rate * time; %mM = mmol/l
% 
% glucose_uptake_flux = glucose_uptake_rate * volume_of_cell; %mmol/cell/hr
% uptaken_glucose_amount_by_flux = glucose_uptake_flux * time * integrated_cell_count; %mmol
% 
% 
% e_0 = glucose_uptake_flux/glu_con;


