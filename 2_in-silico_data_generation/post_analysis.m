%post analysis
%%load max_WT_vs_KRAS_CRC_media\data.mat

%WT_size = size(WT_Model.sol_lst,2);
M =  readmatrix("./data.csv");

glu = zeros(size(M,1),1);
glt = zeros(size(M,1),1);
lac = zeros(size(M,1),1);
biomass = zeros(size(M,1),1);

for i=1:size(M,1)
    glu(i) = M(i,1); %glu
    glt(i) = M(i,2); %glt
    lac(i) = M(i,3); %lac
    biomass(i) = M(i,4); %biomass

end
figure(3)
plot(lac,biomass)
xlabel('lactate secretion rate (mM/hr)')
ylabel('biomass creation rate (1/hr)')

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