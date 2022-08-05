%post analysis
load max_WT_vs_KRAS_CRC_media\data.mat

WT_size = size(WT_Model.sol_lst,2);

glu = zeros(WT_size,1);
glt = zeros(WT_size,1);
lac = zeros(WT_size,1);
biomass = zeros(WT_size,1);

for i=1:WT_size
    glu(i) = WT_Model.sol_lst{1,i}.full(69); %glu
    glt(i) = WT_Model.sol_lst{1,i}.full(70); %glt
    lac(i) = WT_Model.sol_lst{1,i}.full(68); %lac
    biomass(i) = WT_Model.sol_lst{1,i}.full(73); %biomass

end

figure(1)
[X1,Y1] = meshgrid(min(glu):0.0001:max(glu),min(glt):0.0001:max(glt));
Z1=griddata(glu,glt,biomass,X1,Y1);
contourf(X1,Y1,Z1)
xlabel('glucose uptake rate (mM/hr)')
ylabel('glutamine uptake rate (mM/hr)')
title('Biomass creation rate (1/hr)')

figure(2)
[X2,Y2] = meshgrid(min(glu):0.0001:max(glu),min(glt):0.0001:max(glt));
Z2=griddata(glu,glt,lac,X2,Y2);
contourf(X2,Y2,Z2)
xlabel('glucose uptake rate (mM/hr)')
ylabel('glutamine uptake rate (mM/hr)')
title('Lactate creation rate (mM/hr)')