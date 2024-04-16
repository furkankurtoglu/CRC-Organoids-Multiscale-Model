clear
cd output\

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
    growth_rate = log(day_three(k,2)/day_zero(k,2))/72;
    growth_rates(k,:) = [day_three(k,1) growth_rate];
end

insilico_Mean_Growth_rate = mean(growth_rates(:,2));
insilico_Standard_Deviation = std(growth_rates(:,2));

%%
Experimental_Mean_Growth_Rate = 0.01862706034333242;
Experimental_Standard_Deviation = 0.006490257860865904;
figure(1)
hold on
errorbar(1.25,Experimental_Mean_Growth_Rate,Experimental_Standard_Deviation,'r>')
errorbar(1.75,insilico_Mean_Growth_rate,insilico_Standard_Deviation,'k<')
hold off
xlabel('Data Type')
ylabel('Growth Rate (1/hr)')
ylim([0.02 0.032])
title('Growth Rates')
set(gca, 'XTickMode', 'manual', 'XTick', 1.25:0.5:2,  ...
    'XTickLabelMode', 'manual', 'XTickLabel', {'Experimental Growth Rate', 'In silico Growth Rate'},...
     'XLim',[1,2])