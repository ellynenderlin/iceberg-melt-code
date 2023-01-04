%create composite plots of icebersg from different scenes but the same date
clearvars; close all;
region_name = 'Thwaites';
root_dir = '/Users/ellynenderlin/Research/NSF_Antarctic-icebergs/iceberg-melt/'; 
DEM1_date = '20200123'; DEM2_date = '20201007';

%compile the data from the different image footprints
cd([root_dir,region_name]);
% date_dirs = dir([DEM1_date,'*-',DEM2_date,'*']);
date_txts = dir(['*',DEM1_date,'*-',DEM2_date,'*.csv']); dateref = [];
dt = []; xo = []; yo = []; zo = []; Vo = []; xf = []; yf = []; zf = []; Vf = []; po = []; pf = []; coreg_zo = []; coreg_zf = [];
dz = []; dz_sigma = []; dVdt = []; dVdt_uncert = []; draft = []; draft_uncert = []; Asurf = []; Asurf_uncert = []; Asub = []; Asub_uncert = [];
for i = 1:length(date_txts)
    M=readtable(date_txts(i).name); M = table2array(M);
    dt =  [dt; M(:,1)];
    xo =  [xo; M(:,2)]; yo =  [yo; M(:,3)]; zo =  [zo; M(:,4)]; po = [po; M(:,5)]; Vo = [Vo; M(:,6)];
    xf =  [xf; M(:,7)]; yf =  [yf; M(:,8)]; zf =  [zf; M(:,9)]; pf = [pf; M(:,10)]; Vf = [Vf; M(:,11)];
    coreg_zo = [coreg_zo; M(:,12)]; coreg_zf = [coreg_zf; M(:,13)];
    dz = [dz; M(:,14)]; dz_sigma = [dz_sigma; M(:,15)];
    dVdt = [dVdt; M(:,16)]; dVdt_uncert = [dVdt_uncert; M(:,17)];
    draft = [draft; M(:,18)]; draft_uncert = [draft_uncert; M(:,19)];
    Asurf = [Asurf; M(:,20)]; Asurf_uncert = [Asurf_uncert; M(:,21)];
    Asub = [Asub; M(:,22)]; Asub_uncert = [Asub_uncert; M(:,23)];
    dateref = [dateref; repmat(i,length(M(:,1)),1)]; %add a reference to the time-stamped file in case I need to go back and look at results from individual DEM pairs
    clear M;
end
%add a reference number for each iceberg
for i = 1:length(xo)
   if i<10
       berg_ref(i,:) = char(['0',num2str(i)]);
   else
       berg_ref(i,:) = num2str(i);
   end
end
m = dVdt./Asub; H = draft.*(1027./nanmean([po pf],2));
nameref = strfind(date_txts(1).name,'_');
region_abbrev = date_txts(1).name(1:nameref(1)-1);

%plot to figure-out which icebergs need to be re-run
figure; set(gcf,'position',[100 500 1500 600]); date_cmap = colormap(jet(length(date_txts)));
subplot(1,3,1);
for i = 1:length(m)
    pl(dateref(i)) = plot(Asub(i),dVdt(i),'ok','markersize',24,'markerfacecolor','w','markeredgecolor',date_cmap(dateref(i),:)); hold on;
    text(double(Asub(i))-0.03e5,double(dVdt(i)),berg_ref(i,:),'color',date_cmap(dateref(i),:));
end
set(gca,'fontsize',20); xlabel('Submerged area (m^2)','fontsize',20); ylabel('Meltwater flux (m^3/d)','fontsize',20);
grid on; leg = legend(pl,num2str([1:length(date_txts)]')); set(leg,'location','northwest','fontsize',16);
subplot(1,3,2);
for i = 1:length(m)
    plot(H(i),m(i),'ok','markersize',24,'markerfacecolor','w','markeredgecolor',date_cmap(dateref(i),:)); hold on;
    text(double(H(i))-2,double(m(i)),berg_ref(i,:),'color',date_cmap(dateref(i),:));
end
set(gca,'fontsize',20); xlabel('Average iceberg thickness (m)','fontsize',20); ylabel('Melt rate (m/d)','fontsize',20);
grid on;
subplot(1,3,3);
for i = 1:length(m)
    plot(coreg_zo(i)-coreg_zf(i),dVdt(i),'ok','markersize',24,'markerfacecolor','w','markeredgecolor',date_cmap(dateref(i),:)); hold on;
    text(double(coreg_zo(i)-coreg_zf(i)),double(dVdt(i)),berg_ref(i,:),'color',date_cmap(dateref(i),:));
end
set(gca,'fontsize',20); xlabel('\Delta sea-level adjustment (m)','fontsize',20); ylabel('Meltwater flux (m^3/d)','fontsize',20);
grid on;
saveas(gcf,[region_abbrev,'_',DEM1_date,'-',DEM2_date,'_iceberg-coregistration-check_subplots.png'],'png');

%save merged file & delete parent files
T=table(dt,xo,yo,zo,po,Vo,xf,yf,zf,pf,Vf,coreg_zo,coreg_zf,dz,dz_sigma,dVdt,dVdt_uncert,draft,draft_uncert,Asurf,Asurf_uncert,Asub,Asub_uncert);
column_names = ["TimeSeparation", "X_i", "Y_i", "MedianZ_i",...
    "Density_i", "Volume_i", "X_f", "Y_f", "MedianZ_f",...
    "Density_f", "Volume_f", "VerticalAdjustment_i", "VerticalAdjustment_f",...
    "ElevationChange_mean", "ElevationChange_stdev", "VolumeChangeRate", "VolumeChangeRate_uncert",...
    "MedianDraft_mean", "MedianDraft_range", "SurfaceArea_mean", "SurfaceArea_range",...
    "SubmergedArea_mean","SubmergedArea_uncert"];
column_units = ["days", "m", "m", "m",...
    "kg/m^3", "m^3", "m", "m", "m",...
    "kg/m^3", "m^3", "m", "m,"...
    "m", "m", "m^3/d", "m^3/d",...
    "m", "m", "m^2", "m^2",...
    "m^2", "m^2"];
T.Properties.VariableNames = column_names; T.Properties.VariableUnits = column_units;
writetable(T,[root_dir,region_name,'/',region_abbrev,'_',date_txts(1).name(nameref(1)+1:nameref(2)-2),'X_iceberg_meltinfo.csv']);
disp('Text file written');
for i = 1:length(date_txts)
    delete_file = ['delete ',date_txts(i).name]; eval(delete_file);
end


