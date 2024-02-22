close all;
clear all;
clc;
rng(0) % same set of random integers each time simulation is ran
%initCobraToolbox(false);
changeCobraSolver('glpk');

all_reaction_fluxes_percentage=cell(1,74);
tic % record time elapsed for each reaction knockdown

% Simulating all knockdowns
% for jj=1:74
%     disp(['knockout number: ', num2str(jj)])
%     all_reaction_fluxes_percentage{1,jj} = will_code(jj);
%     save('single_knockdown_data.mat')
%     toc
% end

% Simulating just one knockdown 
    all_reaction_fluxes_percentage{1,1} = will_code(2);
    save('testttt.mat')
    toc
