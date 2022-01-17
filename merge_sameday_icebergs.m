%create composite plots of icebersg from different scenes but the same date
clearvars; close all;
region_name = 'Thwaites';
root_dir = '/Users/ellynenderlin/Research/NSF_Antarctic-icebergs/iceberg-melt/'; 
DEM1_date = '20191011'; DEM2_date = '20200123';

%compile the data from the different image footprints
cd_glacier_dir = ['cd ',root_dir,region_name]; eval(cd_glacier_dir);
% date_dirs = dir([DEM1_date,'*-',DEM2_date,'*']);
date_txts = dir(['*',DEM1_date,'*-',DEM2_date,'*.txt']); dateref = [];
dt = []; xo = []; yo = []; zo = []; Vo = []; xf = []; yf = []; zf = []; Vf = []; po = []; pf = []; coreg_zo = []; coreg_zf = [];
dz = []; dz_sigma = []; dVdt = []; dVdt_uncert = []; draft = []; draft_uncert = []; Asurf = []; Asurf_uncert = []; Asub = []; Asub_uncert = [];
for i = 1:length(date_txts)
    load_data = ['M=dlmread(''',date_txts(i).name,''');']; eval(load_data);
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
column_vals = [dt xo yo zo po Vo xf yf zf pf Vf coreg_zo coreg_zf dz dz_sigma dVdt dVdt_uncert draft draft_uncert Asurf Asurf_uncert Asub Asub_uncert dateref];
%merged file will be obvious because second date will have its seconds as X
dlmwrite([region_abbrev,'_',date_txts(1).name(nameref(1)+1:nameref(2)-2),'X_iceberg_meltinfo.txt'],column_vals,'delimiter','\t');
for i = 1:length(date_txts)
    delete_file = ['delete ',date_txts(i).name]; eval(delete_file);
end


%% NO LONGER USE CODE BELOW TO CHECK MERGED DATA... REMOVE
%automatically "fix" melt rate estimates with bad local sea level adjustments
disp('Automatically adjusting fluxes & melt rates for icebergs with clearly bad sea level estimates');
median_sldz = nanmedian(coreg_zo-coreg_zf); median_slzo = nanmedian(coreg_zo); median_slzf = nanmedian(coreg_zf);
rho_i = 900; rho_i_err = 20; %kg m^-3 
rho_sw = 1026;  rho_sw_err = 2; %kg m^-3
for i = 1:length(coreg_zo)
    if (coreg_zo(i)-coreg_zf(i)) > (median_sldz + 1) || (coreg_zo(i)-coreg_zf(i)) < (median_sldz - 1)
        %adjust the volume change estimate so that it is calculated
        %assuming the median sea level correction is accurate locally
        rho_f = (po(i)+pf(i))/2;
        dZ_mean = (rho_sw/(rho_sw-rho_f))*(SL(i).mean.dz-((coreg_zo(i)-coreg_zf(i))-median_sldz));
        dV_mean = SL(i).mean.SA*(dZ_mean+SL(i).SMB+SL(i).creep_dz);
        
        %convert volume change to flux & melt rate
        to = SL(i).initial.time; tf = SL(i).final.time;
        if mod(str2num(to(1:4)),4)==0; doyso=366; modayso = [31 29 31 30 31 30 31 31 30 31 30 31]; else doyso=365; modayso = [31 28 31 30 31 30 31 31 30 31 30 31]; end
        if mod(str2num(tf(1:4)),4)==0; doysf=366; modaysf = [31 29 31 30 31 30 31 31 30 31 30 31]; else doysf=365; modaysf = [31 28 31 30 31 30 31 31 30 31 30 31]; end
        doyo = sum(modayso(1:str2num(to(5:6))))-31+str2num(to(7:8)); doyf = sum(modaysf(1:str2num(tf(5:6))))-31+str2num(tf(7:8));
        if str2num(tf(1:4)) == str2num(to(1:4))
            ddays = doyf-doyo+1;
        elseif str2num(tf(1:4)) - str2num(to(1:4)) == 1
            ddays = doyf + (doyso-doyo)+1;
        else
            years = str2num(to(1:4)):1:str2num(tf(1:4));
            for k = 1:length(years)
                if mod(years(k),4)==0
                    doys(k)=366;
                else
                    doys(k) = 365;
                end
            end
            ddays = doyf + sum(doys(2:end-1)) + (doyso-doyo)+1;
        end
        hrs_o = ((str2num(to(13:14))/(60*60*24))+(str2num(to(11:12))/(60*24))+(str2num(to(9:10))/24));
        hrs_f = ((str2num(tf(13:14))/(60*60*24))+(str2num(tf(11:12))/(60*24))+(str2num(tf(9:10))/24));
        dhrs = hrs_f - hrs_o;
        dt = ddays + dhrs;
        dVdt_mean = dV_mean/dt;
        dHdt_mean = dVdt_mean/SL(i).mean.TA;
        
        %replace in structure
        SL(i).mean.dVdt = dVdt_mean; SL(i).mean.dHdt = dHdt_mean; 
        SL(i).initial.coreg_z = median_slzo; SL(i).final.coreg_z = median_slzf;
    end
end

%replot
dVdt = []; Asub = []; H = []; m = []; coreg_zo = []; coreg_zf = []; berg_ref =[];
for i = 1:length(SL)
    if SL(i).mean.dVdt > 0 && ~isempty(SL(i).mean.TA)
        dVdt = [dVdt SL(i).mean.dVdt];
        Asub = [Asub SL(i).mean.TA];
        H = [H SL(i).mean.H];
        m = [m SL(i).mean.dHdt];
        coreg_zo = [coreg_zo SL(i).initial.coreg_z]; coreg_zf = [coreg_zf SL(i).final.coreg_z];
        berg_ref = [berg_ref; SL(i).name(end-1:end)];
    end
end
figure; set(gcf,'position',[100 100 1200 600]);
subplot(1,2,1);
for i = 1:length(m)
    pl(dateref(i)) = plot(Asub(i),dVdt(i),'ok','markersize',24,'markerfacecolor','w','markeredgecolor',date_cmap(dateref(i),:)); hold on;
    text(double(Asub(i))-0.03e5,double(dVdt(i)),berg_ref(i,:),'color',date_cmap(dateref(i),:));
end
set(gca,'fontsize',20); xlabel('Submerged area (m^2)','fontsize',20); ylabel('Meltwater flux (m^3/d)','fontsize',20);
grid on; leg = legend(pl,dates_string); set(leg,'location','northwest','fontsize',16);
subplot(1,2,2);
for i = 1:length(m)
    plot(H(i),m(i),'ok','markersize',24,'markerfacecolor','w','markeredgecolor',date_cmap(dateref(i),:)); hold on;
    text(double(H(i))-2,double(m(i)),berg_ref(i,:),'color',date_cmap(dateref(i),:));
end
set(gca,'fontsize',20); xlabel('Average iceberg thickness (m)','fontsize',20); ylabel('Melt rate (m/d)','fontsize',20);
grid on;
saveas(gcf,[region_abbrev,'_',DEM1_date,'-',DEM2_date,'_iceberg-melt_subplots.png'],'png');
disp('Iceberg meltwater flux should increase linearly with submerged area');
disp('Iceberg melt rates should increase with thickness');

%call-out the clearly bad icebergs
for i = 1:length(SL)
    if SL(i).mean.dVdt < 0
        disp(['Recalculate elevation change for iceberg #',num2str(i)]);
    end
end        

disp(' ');
disp('If any icebergs stand out as outliers, redo them following the typical steps');
disp(' ');