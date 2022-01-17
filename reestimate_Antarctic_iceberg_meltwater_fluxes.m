%re-extract iceberg meltwater fluxes & melt rates for Antarctic icebergs
%that were pulled prior to updating the interpolation approach
clear all; close all; drawnow;
region_name = 'LarsenB'; region_abbrev = 'LB';
root_dir = '/Users/mariamadryak/Desktop/Antarctic_icebergs/'; 
DEM1_time = '20160117124920'; DEM2_time = '20170107162000';
cd_root_dir = ['cd ',root_dir]; eval(cd_root_dir);
rho_sw = 1026;  rho_sw_err = 2; %kg m^-3
rho_i = 917; rho_i_err = 10; %kg m^-3 

%get DEM names
DEM1_date = DEM1_time(1:8); DEM2_date = DEM2_time(1:8);
glacier_dir = [root_dir,region_name,'/',num2str(DEM1_date),'-',num2str(DEM2_date),'/'];
DEMs = dir([glacier_dir,'/*DEM.mat']);
for i = 1:size(DEMs,1)
    if strmatch(num2str(DEMs(i).name(end-21:end-8)),num2str(DEM1_time))
        DEM1_name = DEMs(i).name;
        disp(['DEM1 = ',num2str(DEM1_name)]);
    elseif strmatch(num2str(DEMs(i).name(end-21:end-8)),num2str(DEM2_time))
        DEM2_name = DEMs(i).name;
        disp(['DEM2 = ',num2str(DEM2_name)]);
    end
end

%pull date/time information
to = DEM1_time; tf = DEM2_time;
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

%load the regional elevation change dataset
cd_to_DEM_offset_data = ['cd ',glacier_dir,'/DEM_offset_data']; eval(cd_to_DEM_offset_data);
load_fjord_data = ['load fjord_offset_',num2str(DEM1_time),'-',num2str(DEM2_time),'.mat']; eval(load_fjord_data);
fjordz_o = fjord.dem1.z(~isnan(fjord.dem1.z)); fjordz_o(fjordz_o>100) = NaN; fjordz_o(fjordz_o<-100) = NaN; 
fjordz_medo = nanmedian(fjord.dem1.z(~isnan(fjord.dem1.z))); fjordz_mado = mad(fjord.dem1.z(~isnan(fjord.dem1.z)),1);
fjordz_o(fjordz_o<fjordz_medo-3*(1.4826*fjordz_mado) | fjordz_o>fjordz_medo+3*(1.4826*fjordz_mado)) = NaN;
fjordz_f = fjord.dem2.z(~isnan(fjord.dem2.z)); fjordz_f(fjordz_f>100) = NaN; fjordz_f(fjordz_f<-100) = NaN; 
fjordz_medf = nanmedian(fjord.dem2.z(~isnan(fjord.dem2.z))); fjordz_madf = mad(fjord.dem2.z(~isnan(fjord.dem2.z)),1);
fjordz_f(fjordz_f<fjordz_medf-3*(1.4826*fjordz_madf) | fjordz_f>fjordz_medf+3*(1.4826*fjordz_madf)) = NaN;

%find the files
iceberg_dir = [glacier_dir,'iceberg_data/'];
cd_to_iceberg_data = ['cd ',iceberg_dir]; eval(cd_to_iceberg_data);
berg_files = dir('*dz.mat');
melt_file = dir('*iceberg_melt.mat');
load_meltfile = ['load ',melt_file(1).name]; eval(load_meltfile);

%for each iceberg, interpolate the elevations from the later date to the
%earlier date
for i = 1:length(berg_files)
    disp(['looping through iceberg ',num2str(i), ' of ',num2str(length(berg_files))]);
    load_file = ['load ',berg_files(i).name]; eval(load_file);
    for j = 1:length(SL)
        if strmatch(berg_files(i).name(1:end-7),SL(j).name(end-8:end))
            SLref = j;
        end
    end
    
    %re-extract elevation change estimates for each translation & rotation iteration
    [xgrid,ygrid] = meshgrid(IB(1).xo,IB(1).yo);
    for p = 1:10
        %no elpsd field (which maintained NaNs in data gaps) so use the
        %difference between the mean & median elevations to identify whether
        %the gaps are filled with anomalously low or high elevations
%         gapo_mask = IB(p).zo.elpsd_adjust.map; gapo_mask(gapo_mask~=0) = 1; gapo_mask(isnan(IB(p).zo.elpsd_adjust.map)) = 0;
        gapo_mask = IB(p).zo.reg_adjust.map;
        if nanmean(IB(p).zo.reg_adjust.map(~isnan(IB(p).zo.reg_adjust.map))) < nanmedian(IB(p).zo.reg_adjust.map(~isnan(IB(p).zo.reg_adjust.map)))
            gapo_mask(IB(p).zo.reg_adjust.map==min(IB(p).zo.reg_adjust.map(~isnan(IB(p).zo.reg_adjust.map)))) = 0;
            gapo_mask(gapo_mask~=0) = 1; gapo_mask(isnan(IB(p).zo.reg_adjust.map)) = 0;
        else
            gapo_mask(IB(p).zo.reg_adjust.map==max(IB(p).zo.reg_adjust.map(~isnan(IB(p).zo.reg_adjust.map)))) = 0;
            gapo_mask(gapo_mask~=0) = 1; gapo_mask(isnan(IB(p).zo.reg_adjust.map)) = 0;
        end
        IB(p).zo.gap_mask = gapo_mask; 
        gapf_mask = IB(p).zf.reg_adjust.map;
        if nanmean(IB(p).zf.reg_adjust.map(~isnan(IB(p).zf.reg_adjust.map))) < nanmedian(IB(p).zf.reg_adjust.map(~isnan(IB(p).zf.reg_adjust.map)))
            gapf_mask(IB(p).zf.reg_adjust.map==min(IB(p).zf.reg_adjust.map(~isnan(IB(p).zf.reg_adjust.map)))) = 0;
            gapf_mask(gapf_mask~=0) = 1; gapf_mask(isnan(IB(p).zf.reg_adjust.map)) = 0;
        else
            gapf_mask(IB(p).zf.reg_adjust.map==max(IB(p).zf.reg_adjust.map(~isnan(IB(p).zf.reg_adjust.map)))) = 0;
            gapf_mask(gapf_mask~=0) = 1; gapf_mask(isnan(IB(p).zf.reg_adjust.map)) = 0;
        end
        IB(p).zf.gap_mask = gapf_mask;
        clear gap*_mask;
        IB(p).zo.reg_adjust.map(IB(p).zo.gap_mask==0) = NaN;
        IB(p).zo.local_adjust.map(IB(p).zo.gap_mask==0) = NaN;
        IB(p).zf.reg_adjust.map(IB(p).zf.gap_mask==0) = NaN;
        IB(p).zf.local_adjust.map(IB(p).zf.gap_mask==0) = NaN;
        
        %interpolate the unrotated & untranslated elevations from the later
        %date to the grid for the earlier date
        [x_rot,y_rot] = unrotate_untranslate_iceberg(IB,p);
        zf_rot = griddata(x_rot,y_rot,double(IB(p).zf.local_adjust.map),xgrid,ygrid);
        IB(p).zf.local_adjust.rotmap = zf_rot; clear zf_rot;
        non_nans = ~isnan(IB(p).zo.local_adjust.map) & ~isnan(IB(p).zf.local_adjust.rotmap);
        IB(p).zo.local_adjust.mean = nanmean(IB(p).zo.local_adjust.map(non_nans));
        IB(p).zf.local_adjust.mean = nanmean(IB(p).zf.local_adjust.rotmap(non_nans));
        
        %pull elevation change info
        DEM_diff.local_adjust = IB(p).zo.local_adjust.map - IB(p).zf.local_adjust.rotmap;
        DEM_diff.local_adjust(DEM_diff.local_adjust == 0) = NaN;
        DEM_diff.local_adjust(isnan(IB(p).zo.local_adjust.map)) = NaN;
        DEM_diff.local_adjust(isnan(IB(p).zf.local_adjust.rotmap)) = NaN;
        IB(p).dz.reg_adjust.map = DEM_diff.local_adjust-(IB(p).local_adjust_o-IB(p).local_adjust_f+nanmedian(min(fjordz_o))-nanmedian(min(fjordz_f)))*ones(size(DEM_diff.local_adjust));
        IB(p).dz.local_adjust.map = DEM_diff.local_adjust;
        IB(p).dz.reg_adjust.mean = nanmean(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.mean = nanmean(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        IB(p).dz.reg_adjust.stdev = nanstd(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.stdev = nanstd(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        IB(p).dz.reg_adjust.median = nanmedian(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.median = nanmedian(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        IB(p).dz.reg_adjust.mad = mad(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)),1);
        IB(p).dz.local_adjust.mad = mad(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)),1);
        IB(p).dz.reg_adjust.max = max(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.max = max(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        IB(p).dz.reg_adjust.min = min(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.min = min(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        non_outliers = IB(p).dz.local_adjust.map<IB(p).dz.local_adjust.median+3*(1.4826*IB(p).dz.local_adjust.mad) & IB(p).dz.local_adjust.map>IB(p).dz.local_adjust.median-3*(1.4826*IB(p).dz.local_adjust.mad);
        IB(p).dz.reg_adjust.filtered_mean = nanmean(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_mean = nanmean(IB(p).dz.local_adjust.map(non_outliers));
        IB(p).dz.reg_adjust.filtered_stdev = nanstd(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_stdev = nanstd(IB(p).dz.local_adjust.map(non_outliers));
        IB(p).dz.reg_adjust.filtered_median = nanmedian(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_median = nanmedian(IB(p).dz.local_adjust.map(non_outliers));
        IB(p).dz.reg_adjust.filtered_mad = mad(IB(p).dz.reg_adjust.map(non_outliers & ~isnan(IB(p).dz.reg_adjust.map)),1);
        IB(p).dz.local_adjust.filtered_mad = mad(IB(p).dz.local_adjust.map(non_outliers & ~isnan(IB(p).dz.local_adjust.map)),1);
        IB(p).dz.reg_adjust.filtered_max = max(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_max = max(IB(p).dz.local_adjust.map(non_outliers));
        IB(p).dz.reg_adjust.filtered_min = min(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_min = min(IB(p).dz.local_adjust.map(non_outliers));
        clear DEM_diff;
    end
    
    %plot a map of the average elevation change
    for p = 1:length(IB)
        dzmap(:,:,p) = IB(p).dz.local_adjust.map;
    end
    figure; set(gcf,'position',[50 50 800 500]);
    subplot(1,2,1);
    imagesc(IB(p).xo,IB(p).yo,IB(p).zo.local_adjust.map); axis xy equal; colormap jet; colorbar;
    title([berg_files(i).name(1:end-7),' initial elevation']);
    subplot(1,2,2);
    imagesc(IB(p).xo,IB(p).yo,nanmean(dzmap,3)); axis xy equal; colormap jet; colorbar;
    set(gca,'clim',[0 nanmean(dzmap(~isnan(dzmap)))+(2*nanstd(dzmap(~isnan(dzmap))))]);
    title([berg_files(i).name(1:end-7),' elevation change']);
    drawnow;

    %resave the individual file
    cd_iceberg_dir = ['cd ',glacier_dir,'/iceberg_data']; eval(cd_iceberg_dir);
    save_file = ['save ',berg_files(i).name,' IB']; eval(save_file);
    
    %pull the elevation change estimates & convert to updated melt rate estimates
    vals = [];
    for j = 1:length(IB)
        if IB(j).local_adjust_f ~= 0 %use local sea level-adjusted elevations
            non_NaNs = length(IB(j).dz.local_adjust.map(~isnan(IB(j).dz.local_adjust.map)));
            vals(size(vals,2)+1:size(vals,2)+non_NaNs) = IB(j).dz.local_adjust.map(~isnan(IB(j).dz.local_adjust.map));
        end
    end
    mean_vals = nanmean(vals); std_vals = nanstd(vals);
    median_vals = nanmedian(vals); mad_vals = mad(vals,1);
    if std_vals > 10*mad_vals %first take at filtering HUGE outliers
        outliers = vals>median_vals+3*2.9 | vals<median_vals-3*2.9; %filter using the median & quoted DEM uncertainty (2.9m)
        vals(outliers) = NaN;
    end
    mean_vals = nanmean(vals); std_vals = nanstd(vals);
    median_vals = nanmedian(vals); mad_vals = mad(vals,1);
    outliers = vals>median_vals+4*1.4826*mad_vals | vals<median_vals-4*1.4826*mad_vals; 
    vals(outliers) = NaN;
    SL(SLref).mean.dz = nanmean(vals); SL(SLref).uncert.dz = nanstd(vals);
    
    %convert to thickness change
    if SL(SLref).orientation == 0
        rho_f = nanmean([SL(SLref).initial.density SL(SLref).final.density]);
        rho_f_range = [rho_i SL(SLref).initial_range.density SL(SLref).final_range.density];
        rho_f_err = [min(rho_f_range-rho_f) max(rho_f_range-rho_f)];
    else
        rho_f = nanmean([SL(SLref).initial.density SL(SLref).final.density]);
        rho_f_range = [SL(SLref).initial_range.density SL(SLref).final_range.density];
        rho_f_err = [min(rho_f_range-rho_f) max(rho_f_range-rho_f)];
    end
    rho_conversion = rho_sw/(rho_sw-rho_f); 
    dZ_mean = (rho_sw/(rho_sw-rho_f))*SL(SLref).mean.dz; 
    dZ_stdev = (rho_sw/(rho_sw-rho_f))*SL(SLref).uncert.dz;
    
    %use elevation change data to get updated estimates & errors
    dH_SMBadjust_mean = dZ_mean + SL(SLref).SMB;
    dH_submelt = dH_SMBadjust_mean + SL(SLref).creep_dz;
    for j = 1:10
        dzs(j) = IB(j).dz.local_adjust.mean;
    end
    SL(SLref).uncert.h_user = nanstd(dzs);
    dH_total_err = abs(dH_submelt).*sqrt(repmat(((SL(SLref).uncert.h_SMB^2 + SL(SLref).uncert.h_rand^2 + SL(SLref).uncert.h_user^2)/SL(SLref).mean.dz^2),size(SL(SLref).uncert.hH_conversion)) + (SL(SLref).uncert.hH_conversion./repmat(rho_conversion,size(SL(SLref).uncert.hH_conversion))).^2);
    dV_mean = SL(SLref).mean.SA*dH_submelt;
    dV_stdev = abs(dV_mean).*sqrt((dZ_stdev./dH_submelt)^2);
    dV_total_err = abs(dV_mean).*sqrt((dH_total_err./repmat(dH_submelt,size(dH_total_err))).^2 + repmat((SL(SLref).change.SA/2)./SL(SLref).mean.SA,size(dH_total_err)).^2);
    dVdt_mean = dV_mean/dt;  dVdt_stdev = dV_stdev/dt;
    dVdt_total_err = dV_total_err./dt;
    SL(SLref).mean.dVdt = dVdt_mean; 
    SL(SLref).range.dVdt(1) = dVdt_mean-3*dVdt_stdev; SL(SLref).range.dVdt(2) = dVdt_mean+3*dVdt_stdev; 
    SL(SLref).uncert.dVdt = dVdt_total_err;
    if ~isempty(SL(SLref).mean.TA)
        dHdt_mean = dVdt_mean/SL(SLref).mean.TA;
        dHdt_stdev = abs(dHdt_mean).*sqrt((dVdt_stdev./dVdt_mean)^2+(SL(SLref).uncert.TA./SL(SLref).mean.TA)^2);
        dHdt_total_err = abs(dHdt_mean).*sqrt((dVdt_total_err./repmat(dVdt_mean,size(dVdt_total_err))).^2 + repmat(SL(SLref).uncert.TA./SL(SLref).mean.TA,size(dVdt_total_err)).^2);
        SL(SLref).mean.dHdt = dHdt_mean;
        SL(SLref).range.dHdt(1) = dHdt_mean-3*dHdt_stdev; SL(SLref).range.dHdt(2) = dHdt_mean+3*dHdt_stdev;
        SL(SLref).uncert.dHdt = dHdt_total_err;
    end
    
    %resave
    cd_to_iceberg_data = ['cd ',iceberg_dir]; eval(cd_to_iceberg_data);
    save_file = ['save(''',melt_file(1).name,''',','''SL''',',','''-v7.3''',')']; eval(save_file);
    clear SLref *vals* outliers rho_f* rho_conversion dz* dZ* dH* dV*;
end

%resave to tab-delimited text file
disp('saving final results to a tab-delimited text file');
cd_to_iceberg_data = ['cd ',iceberg_dir]; eval(cd_to_iceberg_data);
clear dt xo yo zo po Vo xf yf zf pf Vf coreg_z* dz* dVdt* draft* Asurf* Asub*;
column_names = ['dt','xo','yo','zo','po','Vo','xf','yf','zf','pf','Vf','coreg_zo','coreg_zf','dz','dz_sigma','dVdt','dVdt_uncert','draft','draft_uncert','Asurf','Asurf_uncert','Asub','Asub_uncert'];
append_ref = 1;
for i = 1:length(SL)
    if SL(i).mean.dVdt > 0 && ~isempty(SL(i).mean.TA)
        dt(append_ref) = sum(SL(i).days);
        xo(append_ref) = nanmean(SL(i).initial.x); yo(append_ref) = nanmean(SL(i).initial.y); zo(append_ref) = SL(i).initial.z_median; Vo(append_ref) = SL(i).initial.V;
        xf(append_ref) = nanmean(SL(i).final.x); yf(append_ref) = nanmean(SL(i).final.y); zf(append_ref) = SL(i).final.z_median; Vf(append_ref) = SL(i).final.V;
        po(append_ref)=SL(i).initial.density; pf(append_ref) = SL(i).final.density;
        coreg_zo(append_ref) = SL(i).initial.coreg_z; coreg_zf(append_ref) = SL(i).final.coreg_z;
        dz(append_ref) = SL(i).mean.dz; dz_sigma(append_ref) = SL(i).uncert.dz;
        dVdt(append_ref) = SL(i).mean.dVdt; dVdt_uncert(append_ref) = max(SL(i).uncert.dVdt);
        draft(append_ref) = SL(i).mean.draft; draft_uncert(append_ref) = SL(i).change.draft/2;
        Asurf(append_ref) = SL(i).mean.SA; Asurf_uncert(append_ref) = SL(i).change.SA/2;
        Asub(append_ref) = SL(i).mean.TA; Asub_uncert(append_ref) = SL(i).change.TA/2;
        append_ref = append_ref+1;
    else
        disp(['Skipping over data for ',num2str(i)]);
    end
end
column_vals = [dt' xo' yo' zo' po' Vo' xf' yf' zf' pf' Vf' coreg_zo' coreg_zf' dz' dz_sigma' dVdt' dVdt_uncert' draft' draft_uncert' Asurf' Asurf_uncert' Asub' Asub_uncert'];
dlmwrite([region_abbrev,'_',DEM1_time(1:8),'-',DEM2_time(1:8),'_iceberg_meltinfo.txt'],column_vals,'delimiter','\t');
disp('Text file written');

%plot the results
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
figure; set(gcf,'position',[100 500 1500 600]);
subplot(1,3,1);
plot(Asub,dVdt,'ok','markersize',24,'markerfacecolor','w'); hold on;
set(gca,'fontsize',20); xlabel('Submerged area (m^2)','fontsize',20); ylabel('Meltwater flux (m^3/d)','fontsize',20);
for i = 1:length(m)
    text(double(Asub(i))-0.03e5,double(dVdt(i)),berg_ref(i,:))
end
grid on;
subplot(1,3,2);
plot(H,m,'ok','markersize',24,'markerfacecolor','w'); hold on;
set(gca,'fontsize',20); xlabel('Average iceberg thickness (m)','fontsize',20); ylabel('Melt rate (m/d)','fontsize',20);
for i = 1:length(m)
    text(double(H(i))-2,double(m(i)),berg_ref(i,:))
end
grid on;
subplot(1,3,3);
plot(coreg_zo-coreg_zf,dVdt,'ok','markersize',24,'markerfacecolor','w'); hold on;
set(gca,'fontsize',20); xlabel('\Delta sea-level adjustment (m)','fontsize',20); ylabel('Meltwater flux (m^3/d)','fontsize',20);
for i = 1:length(m)
    text(double(coreg_zo(i)-coreg_zf(i)),double(dVdt(i)),berg_ref(i,:))
end
grid on;


