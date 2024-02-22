clear
cd output_spheroid_no_int__NOT_same_size__20\

load output00000000_cells.mat


indices = unique(cells(107,:));
day_zero = zeros(length(indices),2);
day_three = zeros(length(indices),2);


for i = 1:length(indices)
    cell_count = length(find(cells(107,:) == i));
    day_zero(i,1) = indices(i);
    day_zero(i,2) = cell_count;
end

day_zero(end,:) = [];


load output00000012_cells.mat
for i = 1:length(indices)
    cell_count = length(find(cells(107,:) == i));
    day_three(i,1) = indices(i);
    day_three(i,2) = cell_count;
end
day_three(end,:) = [];

growth_rates = zeros(size(day_three,2),2);

for k = 1:size(day_three,1)
    growth_rate =   
    growth_rates(k,:) = [day_three(k,1) growth_rate];
end

insilico_Mean_Growth_rate = mean(growth_rates(:,2));
insilico_Standard_Deviation = std(growth_rates(:,2));
insilico_Median = median(growth_rates(:,2));

%%
Experimental_Growth_Rates = [0.02713954 0.030238357 0.02678475 0.024873669 0.024652633 0.029269302];
Experimental_Mean_Growth_Rate = 0.017618024156208213;
Experimental_Standard_Deviation = 0.008428285380733428;
Experimental_Median_Growth_Rate = 0.026962145;
ExperimentalQR = prctile(Experimental_Growth_Rates,[25 75]);
Experimental_LowerDistance = Experimental_Median_Growth_Rate - ExperimentalQR(1);
Experimental_HigherDistance = ExperimentalQR(2) - Experimental_Median_Growth_Rate;

figure(1)
hold on
errorbar(1.25,Experimental_Mean_Growth_Rate,Experimental_Standard_Deviation,'k>','LineWidth',2)
errorbar(1.75,insilico_Mean_Growth_rate,insilico_Standard_Deviation,'r<','LineWidth',2)
hold off
ylabel('Growth Rate (1/hr)')
ylim([0.02 0.032])
title('Spheroid Growth Rates')
set(gca, 'XTickMode', 'manual', 'XTick', 1.25:0.5:2,  ...
    'XTickLabelMode', 'manual', 'XTickLabel', {'Experimental', 'Simulation'},...
     'XLim',[1.2,1.80])
