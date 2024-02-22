function Met_fluxes_percents = will_code(j)
%run_reaction_fcn runs the reaction for a given j (knockdown reaction) and
%returns Met_fluxes_percents
% input(s): j(int) - knockdown reaction number
% output(s): Met_fluxes_percents(cell) - flux of the knockedout reaction
save_data = 'Y';




Met_numpoints=8; %Number of evenly spaced points of enzyme knockout percentage from 0%-100%, eg. if numpoints=4, 0% 33% 66% 100%
Met_fluxes_percents=cell(1,Met_numpoints);
Met_bio_percents=cell(1,Met_numpoints);
Met_percents_values=zeros(1,Met_numpoints); %Matrix that saves percentages of enzyme knockouts

for y=1:Met_numpoints
    percent = round(((y*(1/(Met_numpoints-1))) - (1/(Met_numpoints-1))),1); %percent in decimal form of enzyme knockout, saved in matrix
    Met_percents_values(1,y)=percent;


%folder_path_prefix = '';  % specify the location of "data.mat" (flux distributions of the original phenotype)
% folder_path_lst = {[folder_path_prefix 'WT_vs_KRAS_CRC_media/'], [folder_path_prefix 'WT_vs_KRAS_CAF_media/'], ...
% [folder_path_prefix 'WT_CRC_media_vs_CAF_media/'], [folder_path_prefix 'KRAS_CRC_media_vs_CAF_media/']};  % these four folders should have been created. 

folder_path_lst = {['baselineData/WT_vs_KRAS_CRC_media/'], ['baselineData/WT_vs_KRAS_CAF_media/'], ...
    ['baselineData/WT_CRC_media_vs_CAF_media/'], ['baselineData/KRAS_CRC_media_vs_CAF_media/']};  


%Getting matrix of fluxes for each condition
for z=1:2
    if z==1
        %load('KRAS_CRC_media_vs_CAF_media/data.mat');
        load('baselineData/KRAS_CRC_media_vs_CAF_media/data.mat');
        KRAS_CRC=data(:,1:2:end);
        KRAS_CCM=data(:,2:2:end);
    end
    if z==2
        %load('WT_CRC_media_vs_CAF_media/data.mat');
        load('baselineData/WT_CRC_media_vs_CAF_media/data.mat');
        WT_CRC=data(:,1:2:end);
        WT_CCM=data(:,2:2:end);
    end
end

b_nonzero_indices_mat = zeros(4, 100, 74, 89);  % 4 conditions, 100 sets of mass balance constraints, 74 reactions, 89 metabolites
knockdown_biomass_mat = cell(4, 100, 1);
%for i=1:4    % i=1, WT(CRC-only media); i=2, KRAS(CRC-CAF media); i=3, WT(CRC-CAF media); i=4, KRAS(CRC-only media), this is the correct order
i=1;
data_path = [folder_path_lst{i} 'data.mat'];
load(data_path);

for a=1:Met_numpoints
    for b = 1:Met_numpoints
        for k=1:100
            glucose_percent = (1-(a-1)/Met_numpoints);

            glutamine_level = 0.003*(1-(b-1)/Met_numpoints);
            %disp(glucose_percent)
            if (i == 2) || (i == 3)
                model = KRAS_Model.model_lst{1, k};
                model.b(82) = KRAS_Model.v_lst(73, k);
            else
                model = WT_Model.model_lst{1, k};
                model.b(82) = WT_Model.v_lst(73, k);
            end
            model.b(72) = 0;
            %4 if statements that give fractional enzyme knockouts depending on
            %condition being used
            if i==1
                if mean(WT_CRC(:,j)) > 0
                    tmp_model = changeRxnBounds(model, model.rxns{j}, (mean(WT_CRC(:,j)))*(1-percent), 'u');
                    if glucose_percent > (1-percent)
                        glucose_level = 0.223*(1-percent);
                        tmp_model = changeRxnBounds(tmp_model,'GLUT',glucose_level,'b');
                    else
                        glucose_level = 0.223*glucose_percent;
                        tmp_model = changeRxnBounds(tmp_model,'GLUT',glucose_level,'b');
                    end
                    tmp_model = changeRxnBounds(tmp_model,'ASCT2',glutamine_level,'b');
%                    tmp_model = changeRxnBounds(tmp_model, tmp_model.rxns{j}, 0.0, 'l');
                    tmp_sol = optimizeCbModel(tmp_model, 'max');
            end
        %     if i==2
        %         tmp_model = changeRxnBounds(model, model.rxns{j}, (mean(KRAS_CCM(:,j)))*(1-percent), 'u');
        %         tmp_model = changeRxnBounds(tmp_model, tmp_model.rxns{j}, 0.0, 'l');
        %         tmp_sol = optimizeCbModel(tmp_model, 'max');
        %     end
        %     if i==3
        %         tmp_model = changeRxnBounds(model, model.rxns{j}, (mean(WT_CCM(:,j)))*(1-percent), 'u');
        %         tmp_model = changeRxnBounds(tmp_model, tmp_model.rxns{j}, 0.0, 'l');
        %         tmp_sol = optimizeCbModel(tmp_model, 'max');
        %     end
        %     if i==4
        %         tmp_model = changeRxnBounds(model, model.rxns{j}, (mean(KRAS_CRC(:,j)))*(1-percent), 'u');
        %         tmp_model = changeRxnBounds(tmp_model, tmp_model.rxns{j}, 0.0, 'l');
        %         tmp_sol = optimizeCbModel(tmp_model, 'max');
        %     end
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
                        knockdown_biomass_mat{i, k, j} = NaN;
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
                               knockdown_biomass_mat{i, k, j} = solutionDel.v;
                           else
                               knockdown_biomass_mat{i, k, j} = 0;
                           end
                        else
                           v = zeros(size(tmp_model.rxns));
                           r = zeros(size(tmp_model.mets));
                           p = zeros(size(tmp_model.rxns));
                           q = zeros(size(tmp_model.rxns));
                           knockdown_biomass_mat{i, k, j} = NaN;
                        end
                    end
                else
                    if tmp_sol.f > 0
                        knockdown_biomass_mat{i, k, j} = tmp_sol.v;
                    else
                        knockdown_biomass_mat{i, k, j} = 0;
                    end
                end
                
                %tmp_sol.f
            if save_data == "Y"    
                features_to_save = [glucose_level,glutamine_level,percent,tmp_sol.f];
                %%writematrix(sol_vector,'./data.csv',Delimiter=',',WriteMode='append');
                %writematrix(sol_vector,'../DNN_model_generation/WT_in_silico_data_even_NO_INTRACELLULAR.csv',Delimiter=',',WriteMode='append');
                writematrix(features_to_save,'../3_DNN_model_generation/HK_Gene_Knockout_trials_Glc_Gln_added.csv',Delimiter=',',WriteMode='append')
            end
        end
    end
end
% save("enzyme_knockout_result.mat", "knockdown_biomass_mat", "b_nonzero_indices_mat");

knockdown_biomass_mat_Met=knockdown_biomass_mat; %Rename matrix

%Save Biomass production and put into matrix
A=knockdown_biomass_mat_Met(:,:,j);
B=zeros(4,100);
%for l=1:4
l=1;
for m=1:100
    if A{l,m}==0
        B(l,m)=0;
    else
        B(l,m)=A{l,m}(73);
    end
end

Met_bio_percents{1,y}=B';

%Save fluxes and put into matrix
C=zeros(4,100,74);
m=1;
for m=1:100
    for n=1:74
        if A{l,m}==0
            C(l,m,n)=0;
        else
            C(l,m,n)=A{l,m}(n);
        end
    end
    end

Met_fluxes_percents{1,y}=squeeze(mean(C,2));
end

end

