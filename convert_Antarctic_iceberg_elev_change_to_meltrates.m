function [SL] = convert_Antarctic_iceberg_elev_change_to_meltrates(DEM1,DEM2,IM1,IM2,berg_numbers,region_name,region_abbrev,dir_output,dir_code,dir_iceberg,dir_bedrock)
% Function to convert iceberg elevation change to melt rates in Antarctica
% Ellyn Enderlin & Rainey Aberle, Fall 2021
%
% INPUTS:   DEM1            structure variable containing earlier DEM info
%           DEM2            structure variable containing later DEM info
%           IM1             structure variable containing earlier
%                           orthoimage info
%           IM2             structure variable containing later orthoimage
%                           info
%           berg_numbers    number of iceberg to detect elevation change
%           dir_output      directory where all output files will be placed
%           region_abbrev   abbrevation of region used in file names
%
% OUTPUTS:  SL 

%set densities: assume steady-state profiles (like the Ligtenberg et al. (2014) paper
%describing these data), as supported by fairly constant AIS climate over the past ~40 years
rho_sw = 1026;  rho_sw_err = 2; %kg m^-3
rho_i = 917; rho_i_err = 10; %kg m^-3 

coord_files = dir([dir_iceberg,'*PScoords.txt']);
for j = 1:length(coord_files)
    coords=readmatrix(coord_files(j).name);
    PSy_early(j)=coords(1);  PSx_early(j)=coords(2); PSy_late(j)=coords(3); PSx_late(j)=coords(4);
    clear coords;
end
berg_x = nanmean([PSx_early PSx_late]); berg_y = nanmean([PSy_early PSy_late]);
[berg_lon,berg_lat] = ps2wgs(berg_x,berg_y,'StandardParallel',-71,'StandardMeridian',0);

%read the RACMO files & find nearest pixel
cd([dir_code,'general/RACMO2.3_Antarctica/']);
prompt = 'Are you looking at icebergs along the Antarctic Peninsula (y/n)?';
str = input(prompt,'s');
if strmatch(str,'y')==1
    AP=1;
    runoff_lat = ncread('RACMO2.3p2_XPEN055_runoff_daily_2011_2016.nc','lat');
    runoff_lon = ncread('RACMO2.3p2_XPEN055_runoff_daily_2011_2016.nc','lon');
    runoff = ncread('RACMO2.3p2_XPEN055_runoff_daily_2011_2016.nc','runoff'); runoff(runoff==-9999)=NaN;
    snowmelt = ncread('RACMO2.3p2_XPEN055_snowmelt_daily_2011_2016.nc','snowmelt'); snowmelt(snowmelt==-9999)=NaN;
    airtemp = ncread('RACMO2.3p2_XPEN055_T2m_daily_2011_2016.nc','t2m'); airtemp(airtemp==-9999)=NaN;
    icetemp = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','t2m'); icetemp(icetemp==-9999)=NaN; %assume the ice temp matches the avg long-term temp (for creep estimates)
    smb = [];
    runoff_time = ncread('RACMO2.3p2_XPEN055_runoff_daily_2011_2016.nc','time');
else
    AP=0;
    runoff_lat = ncread('RACMO2.3p2_ANT27_runoff_daily_2011_2016.nc','lat');
    runoff_lon = ncread('RACMO2.3p2_ANT27_runoff_daily_2011_2016.nc','lon');
    runoff = ncread('RACMO2.3p2_ANT27_runoff_daily_2011_2016.nc','runoff'); runoff(runoff==-9999)=NaN;
    snowmelt = [];
    airtemp = ncread('RACMO2.3p2_ANT27_T2m_daily_2011_2016.nc','t2m'); airtemp(airtemp==-9999)=NaN;
    icetemp = ncread('RACMO2.3p2_ANT27_T2m_monthly_1979_2016.nc','t2m'); icetemp(icetemp==-9999)=NaN; %assume the ice temp matches the avg long-term temp (for creep estimates)
    smb = ncread('RACMO2.3p2_ANT27_smb_daily_2011_2016.nc','smb'); smb(smb==-9999)=NaN;
    runoff_time = ncread('RACMO2.3p2_ANT27_runoff_daily_2011_2016.nc','time');
end
runoff_days = runoff_time-runoff_time(1); %days since 20110101
for i = 1:size(runoff_lat,1)
    for j = 1:size(runoff_lat,2)
        [runoff_x(i,j),runoff_y(i,j)] = wgs2ps(runoff_lon(i,j),runoff_lat(i,j),'StandardParallel',-71,'StandardMeridian',0);
    end
end
%calculate the x&y distance between the target of interest (berg lon, berg lat) and each RACMO grid cell
lat_diff = abs(berg_lat - runoff_lat);
lon_diff = abs(berg_lon - runoff_lon);
diff_map = sqrt(lat_diff.^2+lon_diff.^2); diff_map(isnan(squeeze(nanmean(runoff(:,:,1,270:450),4)))) = NaN;%solve for the distance vector using the x&y distances
RACMO_ref = find(diff_map==min(min(diff_map))); %find the minimum distance (reference for your grid cell is output)
[RACMOy RACMOx] = ind2sub(size(squeeze(nanmean(runoff(:,:,1,270:450),4))),RACMO_ref); %convert cell reference to an x- and y-cell index
disp(['RACMO x-reference = ',num2str(RACMOx),' & y-reference = ',num2str(RACMOy)]);

%adjust RACMO reference grid cell if necessary
disp('Adjust coordinates (if necessary) to extract surface melt estimates');
figure; set(gcf,'position',[500 100 900 900]);
runoff_cmap = colormap(jet(10001)); runoff_cmap(1,:) = [1 1 1];
if ~isempty(snowmelt)
    imagesc(max(snowmelt(:,:,270:450),[],3).*86500./1000); colormap(gca,runoff_cmap); hold on; set(gca,'ydir','reverse'); %plot coordinates on average melt from peak summer melt rates in the big 2012 melt season
    disp(['snowmelt at RACMO pixel = ',num2str(max(snowmelt(RACMOx,RACMOy,270:450)).*86500./1000),' mm w.e.']);
else
    imagesc(max(smb(:,:,270:450),[],3).*86500./1000); colormap(gca,runoff_cmap); hold on; set(gca,'ydir','reverse'); %plot coordinates on average melt from peak summer melt rates in the big 2012 melt season
    disp(['surface mass balance at RACMO pixel = ',num2str(max(smb(RACMOx,RACMOy,270:450)).*86500./1000),' mm w.e.']);
end
set(gca,'clim',[0 0.05]); cbar = colorbar; set(get(cbar,'ylabel'),'string','Runoff (m w.e. per day');
plot(RACMOx,RACMOy,'ok','markerfacecolor','k','markeredgecolor','w','markersize',20); hold on;
plot(RACMOx,RACMOy,'xw','markersize',20); hold on;
set(gca,'xlim',[RACMOx-10 RACMOx+10],'ylim',[RACMOy-10 RACMOy+10]);
prompt = 'Modify coordinates if the marker is in a region with no data. Do the coordinates need to be modified (y/n)?';
str = input(prompt,'s');
if strmatch(str,'y')==1
    disp('To change coordinates, identify new coordinates using tick marks and type RACMOx=XX; RACMOy=YY; dbcont ');
    keyboard
end
close all; drawnow;

%extract FAC data from the Ligtenberg model output
cd([dir_code,'general/FDM_Antarctica/']);
%FAC = firn air content (H_observed-FAC = H_i)
if AP==1
    firn_lat = ncread('FDM_FirnAir_XPEN055_1979-2016.nc','lat');
    firn_lon = ncread('FDM_FirnAir_XPEN055_1979-2016.nc','lon');
    FAC = ncread('FDM_FirnAir_XPEN055_1979-2016.nc','FirnAir'); 
    firnu_lat = ncread('FDM_FirnAir-uncertainty_XPEN055_1979-2016.nc','lat'); 
    firnu_lon = ncread('FDM_FirnAir-uncertainty_XPEN055_1979-2016.nc','lon'); 
    FACu = ncread('FDM_FirnAir-uncertainty_XPEN055_1979-2016.nc','Total_unc'); 
    firn_time = ncread('FDM_FirnAir_XPEN055_1979-2016.nc','time');
else
    firn_lat = ncread('FDM_FirnAir_ANT27_1979-2016.nc','lat');
    firn_lon = ncread('FDM_FirnAir_ANT27_1979-2016.nc','lon');
    FAC = ncread('FDM_FirnAir_ANT27_1979-2016.nc','FirnAir'); 
    firnu_lat = ncread('FDM_FirnAir-uncertainty_ANT27_1979-2016.nc','lat'); 
    firnu_lon = ncread('FDM_FirnAir-uncertainty_ANT27_1979-2016.nc','lon'); 
    FACu = ncread('FDM_FirnAir-uncertainty_ANT27_1979-2016.nc','Uncert_total'); 
    firn_time = ncread('FDM_FirnAir_ANT27_1979-2016.nc','time');
end
for i = 1:size(firn_lat,1)
    for j = 1:size(firn_lat,2)
        [firn_x(i,j),firn_y(i,j)] = wgs2ps(firn_lon(i,j),firn_lat(i,j),'StandardParallel',-71,'StandardMeridian',0);
    end
end
clear *diff*;
lat_diff = abs(berg_y*ones(size(firn_y)) - firn_y);
lon_diff = abs(berg_x*ones(size(firn_x)) - firn_x);
diff_map = sqrt(lat_diff.^2+lon_diff.^2);
firn_ref = find(diff_map==min(min(diff_map)));
[firny firnx] = ind2sub(size(nanmean(FAC(:,:,size(FAC,3)-37:size(FAC,3)),3)),firn_ref);
disp(['FAC x-reference = ',num2str(firnx),' & y-reference = ',num2str(firny)]);

%adjust FAC reference grid cell if necessary
disp('Adjust coordinates (if necessary) to extract firn density estimates');
figure; set(gcf,'position',[500 100 900 900]);
imagesc(nanmean(FAC,3)); colormap(jet(10001)); hold on; set(gca,'ydir','reverse'); 
cbar = colorbar; set(get(cbar,'ylabel'),'string','FAC (m) ');
plot(firnx,firny,'ok','markerfacecolor','k','markeredgecolor','w','markersize',20); hold on;
plot(firnx,firny,'xw','markersize',20); hold on;
set(gca,'xlim',[firnx-10 firnx+10],'ylim',[firny-10 firny+10]);
disp(['firn air content estimate = ',num2str(nanmean(FAC(firny,firnx,:),3)),'m']);
prompt = 'Modify coordinates if the marker is in a region with no data. Do the coordinates need to be modified (y/n)?';
str = input(prompt,'s');
if strmatch(str,'y')==1
    disp('To change coordinates, identify new coordinates using tick marks and type firnx=XX; firny=YY; dbcont ');
    keyboard
end
close all; drawnow;

%interpolate the FAC uncertainty map to the FAC grid
for i = 1:size(firnu_lat,1)
    for j = 1:size(firnu_lat,2)
        [firnu_x(i,j),firnu_y(i,j)] = wgs2ps(firnu_lon(i,j),firnu_lat(i,j),'StandardParallel',-71,'StandardMeridian',0);
    end
end
facu_x=reshape(firnu_x,[(size(FACu,1)*size(FACu,2)),1]);
facu_y=reshape(firnu_y,[(size(FACu,1)*size(FACu,2)),1]);
facu=reshape(FACu,[(size(FACu,1)*size(FACu,2)),1]);
F=scatteredInterpolant(facu_x,facu_y,facu);
FACu_interp = F(firn_x,firn_y);

%pull the annual avg FAC & extrapolate density depth profiles to estimate
%uncertainties in iceberg density
firnair.xref=firnx; firnair.yref=firny;
firnair.mean = nanmean(FAC(firny,firnx,:)); firnair.std = nanstd(FAC(firny,firnx,:));
firnair.median = nanmedian(FAC(firny,firnx,:)); firnair.mad = mad(FAC(firny,firnx,:),1);
firnair.uncert = FACu_interp(firny,firnx);
if isnan(firnair.uncert)
    F=scatteredInterpolant(facu_x,facu_y,facu,'nearest');
    FACu_interp = F(firn_x,firn_y);
    firnair.uncert = FACu_interp(firny,firnx);
end
density_lat = ncread('FDM_DepthDensLevels_ANT27_map.nc','lat');
density_lon = ncread('FDM_DepthDensLevels_ANT27_map.nc','lon');
sevenhundred = ncread('FDM_DepthDensLevels_ANT27_map.nc','z700');
sevenfifty = ncread('FDM_DepthDensLevels_ANT27_map.nc','z750');
eighthundred = ncread('FDM_DepthDensLevels_ANT27_map.nc','z800');
eightthirty = ncread('FDM_DepthDensLevels_ANT27_map.nc','z830');
for i = 1:size(density_lat,1)
    for j = 1:size(density_lat,2)
        [density_x(i,j),density_y(i,j)] = wgs2ps(density_lon(i,j),density_lat(i,j),'StandardParallel',-71,'StandardMeridian',0);
    end
end
clear *diff*;
lat_diff = abs(berg_y*ones(size(density_y)) - density_y);
lon_diff = abs(berg_x*ones(size(density_x)) - density_x);
diff_map = sqrt(lat_diff.^2+lon_diff.^2);
density_ref = find(diff_map==min(min(diff_map)));
[densityy,densityx] = ind2sub(size(eightthirty),density_ref);
disp(['Density profile x-reference = ',num2str(densityx),' & y-reference = ',num2str(densityy)]);
disp(['Pore closure depth = ',num2str(eightthirty(densityy,densityx)),'m']);
if isnan(eightthirty(densityy,densityx))
    figure; set(gcf,'position',[500 100 900 900]);
    imagesc(eightthirty); hold on; 
    cmap = colormap(jet(10001)); cmap(1,:) = [1 1 1]; colormap(gca,cmap); set(gca,'ydir','reverse');
    cbar = colorbar; set(get(cbar,'ylabel'),'string','FAC (m) ');
    plot(densityx,densityy,'ok','markerfacecolor','k','markeredgecolor','w','markersize',20); hold on;
    plot(densityx,densityy,'xw','markersize',20); hold on;
    set(gca,'xlim',[densityx-10 densityx+10],'ylim',[densityy-10 densityy+10]);
    disp('To change coordinates, identify new coordinates and type densityx=XX; densityy=YY; dbcont ');
    keyboard
    close all; drawnow;
end
disp(['Pore closure depth = ',num2str(eightthirty(densityy,densityx)),'m']);
density.xref = densityx; density.yref = densityy;
density.sevhun = sevenhundred(density.yref,density.xref);
density.sevfif = sevenfifty(density.yref,density.xref);
density.eighthun = eighthundred(density.yref,density.xref);
density.eightthir = eightthirty(density.yref,density.xref);
density_levels = [700 750 800 830]; density_depths = [density.sevhun density.sevfif density.eighthun density.eightthir];
%fit a curve to the density-depth profile to estimate the depth of the base of the firn column
% [f,gof] = fit(density_levels',density_depths','smoothingspline');
% ex_density = fnxtr(f.p); %extrapolate outside the data domain using a second order polynomial
% density.nineseventeen = fnval(ex_density,917); %evaluate the extrapolation function at the bubble-free ice density
ft = fittype('917-(917-a)*exp(-x/b)');
[f,~] = fit(density_depths',density_levels',ft,'StartPoint',[350 firnair.median]); %empirical, exponential density-depth relation from Schytt (1958) on Cuffey & Paterson p. 19
density.nineseventeen = -f.b*log(-(916.9-917)/(917-f.a)); %find depth where rho=916.9 (goes to infinity at 917)
% density_levels = [density_levels 917]; density_depths = [density_depths density.nineseventeen];

%save the FAC & density data
if ~exist([dir_output,'firn_data/'],'dir')
    mkdir([dir_output,'firn_data/']);
end
save([dir_output,'firn_data/',region_name,'_density_data.mat'],'firnair','density');
close all; drawnow;

%locate the bedrock adjustment file
bedrock_file = dir([dir_bedrock,'bedrock_offset*.mat']);

%calculate the time separation between DEMs in terms of 
%decimal years (ddays) & decimal days (days)
to = DEM1.time; tf = DEM2.time;
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
    
    %calculate the sum of days differently if during a leap year
    if doyo > sum(modayso(1:2))
        ddays = doyf + sum(doys(2:end-1)) + (365-doyo)+1;
    else
        ddays = doyf + sum(doys(2:end-1)) + (366-doyo)+1;
    end
end
hrs_o = ((str2num(to(13:14))/(60*60*24))+(str2num(to(11:12))/(60*24))+(str2num(to(9:10))/24));
hrs_f = ((str2num(tf(13:14))/(60*60*24))+(str2num(tf(11:12))/(60*24))+(str2num(tf(9:10))/24));
dhrs = hrs_f - hrs_o;
dt = ddays + dhrs;
days = ones(1,ceil(dt)); days(1) = 1-hrs_o; days(2:end-1) = 1; days(end) = hrs_f;

%estimate surface melting using RACMO runoff (mm w.e. per day)
%data start in 2011 so date indices are referenced to Jan 1, 2011
clear doys;
yrs = 2011:1:str2num(tf(1:4));
for k = 1:length(yrs)
    if mod(yrs(k),4)==0
        doys(k)=366;
    else
        doys(k) = 365;
    end
end
doy1 = sum(doys(1:str2num(to(1:4))-2011))+doyo;
doy2 = sum(doys(1:str2num(tf(1:4))-2011))+doyf; 
decidayo =  str2num(to(1:4))+(doyo/doyso); decidayf =  str2num(tf(1:4))+(doyf/doysf); 
if decidayf<(2017+(max(runoff_days)-365-366-365-365-365-366)/365)
    melt = squeeze(runoff(RACMOx,RACMOy,doy1:doy2));
else
    if doy1 <= length(runoff)
        melt = cat(1,squeeze(runoff(RACMOx,RACMOy,doy1:length(runoff))),squeeze(runoff(RACMOx,RACMOy,sum(doys(1:2016-2011)):sum(doys(1:2016-2011))+doyf-1))); %temporarily use 2016 melt estimates from the same DOY until RACMO data are updated
    else
        if str2num(tf(1:4)) - str2num(to(1:4)) <= 1
            melt = cat(1,squeeze(runoff(RACMOx,RACMOy,sum(doys(1:2016-2011))+doyo:length(runoff))),squeeze(runoff(RACMOx,RACMOy,sum(doys(1:2016-2011)):sum(doys(1:2016-2011))+doyf-1)));
        else
            melt = cat(1,squeeze(runoff(RACMOx,RACMOy,sum(doys(1:2016-2011))+doyo:length(runoff))),squeeze(runoff(RACMOx,RACMOy,sum(doys(1:2016-2011)):sum(doys(1:2017-2011)))),squeeze(runoff(RACMOx,RACMOy,sum(doys(1:2016-2011)):sum(doys(1:2016-2011))+doyf-1)));
        end
    end
end
if length(days)~=length(melt)
    clear days;
    days = ones(1,length(melt)); days(1) = 1-hrs_o; days(2:end-1) = 1; days(end) = hrs_f;
end
surfmelt = nansum(days'.*melt)/1000; %surface meltwater that runs off (mm w.e. per day)
%estimate the ice temperature as the average long-term air temperature 
%(doesn't account for advection of colder ice and melt/refreezing at the surface and/or submarine interface)
iceberg_avgtemp = nanmean(squeeze(icetemp(RACMOx,RACMOy,:))); % air temp (Kelvin)
% iceberg_stdtemp = nanstd(squeeze(icetemp(RACMOx,RACMOy,:))); % air temp variability (Kelvin)


% %load the saved data if restarting
% cd iceberg_data
% load_saved_data = ['load ',region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat']; eval(load_saved_data);

% load ROI and elevation info for all icebergs
cd(dir_iceberg);
disp('Compiling iceberg polygon information...')
for i = 1:size(berg_numbers,1)
    % loop through DEM data & iceberg outlines
    berg_number = berg_numbers(i).name(1:end-7);
    load_file = ['load ',berg_numbers(i).name]; eval(load_file);
    
    %incorporate data into a structure
    SL(i).name = [region_name,num2str(berg_number)];
    xo = []; yo = []; xf = []; yf = [];
    for j = 1:10
        xo = cat(1,xo,IB(j).vertices.xo); yo = cat(1,yo,IB(j).vertices.yo);
        xf = cat(1,xf,IB(j).vertices.xf); yf = cat(1,yf,IB(j).vertices.yf);   
        non_nans = ~isnan(IB(j).zo.local_adjust.map) & ~isnan(IB(j).zf.local_adjust.rotmap);
        IB(j).zo.local_adjust.map(~non_nans) = NaN; IB(j).zf.local_adjust.rotmap(~non_nans) = NaN;
        zo(:,:,j) = IB(j).zo.local_adjust.map; zf(:,:,j) = IB(j).zf.local_adjust.rotmap;
        clear non_nans;
    end
    SL(i).initial.x = nanmean(xo,1); SL(i).initial.y = nanmean(yo,1);
    SL(i).final.x = nanmean(xf,1); SL(i).final.y = nanmean(yf,1);
%     poly1 = roipoly(DEM1.x,DEM1.y,DEM1.z_elpsd_adjust,SL(i).initial.x,SL(i).initial.y);
%     mask1 = logical(poly1);
%     SL(i).initial.z = single(mask1.*DEM1.z_elpsd_adjust); SL(i).initial.z(SL(i).initial.z == 0) = NaN; SL(i).initial.z(SL(i).initial.z<0) = NaN;
    SL(i).initial.z = nanmean(zo,3);
%     poly2 = roipoly(Y.x,Y.y,Y.z_elpsd_adjust,SL(i).final.x,SL(i).final.y);
%     mask2 = logical(poly2);
%     SL(i).final.z = single(mask2.*Y.z_elpsd_adjust); SL(i).final.z(SL(i).final.z == 0) = NaN; SL(i).final.z(SL(i).final.z<0) = NaN;
    SL(i).final.z = nanmean(zf,3);
    SL(i).initial.coreg_z = IB(1).local_adjust_o; SL(i).final.coreg_z = IB(1).local_adjust_f;
    SL(i).initial.time = to; SL(i).final.time = tf; SL(i).days = days;
    SL(i).SMB = -surfmelt; SL(i).airtemp = iceberg_avgtemp;
    if isfield(IB,'flag')
        SL(i).flag = IB(1).flag;
    else
        SL(i).flag = NaN;
    end
    clear zo zf;
end
disp('...done');
clear Y Z;

% ----------load the initial image & DEM files to map the initial iceberg margins----------
%set image bounds & crop
xo = []; yo = [];
for i = 1:length(SL)
    xo = [xo SL(i).initial.x]; yo = [yo SL(i).initial.y];
end
dy = IM1.y(1)-IM1.y(2);
if dy < 0
    x1 = find(IM1.x >= min(xo)-1500,1,'first'); x2 = find(IM1.x <= max(xo)+1500,1,'last');
    y1 = find(IM1.y >= max(yo)+1500,1,'first'); y2 = find(IM1.y <= min(yo)-1500,1,'last');
    xlims = sort([x1, x2]); ylims = sort([y1, y2]);
    xmin=xlims(1); xmax=xlims(2);
    xmin(isempty(xmin)) = 1; xmax(isempty(xmax)) = length(IM1.x);
    if length(ylims) ==1
        y2(isempty(y2)) = 1; y1(isempty(y1)) = length(IM1.y);
    end
    ylims = sort([y1, y2]);
    ymin=ylims(1); ymax=ylims(2);
else
    x1 = find(IM1.x >= min(xo)-1500,1,'first'); x2 = find(IM1.x <= max(xo)+1500,1,'last');
    y1 = find(IM1.y >= max(yo)+1500,1,'last'); y2 = find(IM1.y <= min(yo)-1500,1,'first');
    xlims = sort([x1, x2]); ylims = sort([y1, y2]);
    xmin=xlims(1); xmax=xlims(2);
    xmin(isempty(xmin)) = 1; xmax(isempty(xmax)) = length(IM1.x);
    if length(ylims)==1 
        y1(isempty(y1)) = 1; y2(isempty(y2)) = length(IM1.y);
    end
    ylims = sort([y1, y2]);
    ymin=ylims(1); ymax=ylims(2);
end
A.x = IM1.x(1,[xlims(1):xlims(2)]); A.y = IM1.y(1,ylims(1):ylims(2));
A.z = double(IM1.z(ylims(1):ylims(2),xlims(1):xlims(2)));
if mean(A.z(~isnan(A.z)),'omitnan') < (1/3)*255 %image is dark
    z_adjust = imadjust(A.z./255,[],[],0.25);
    A.z = z_adjust; clear z_adjust;
end

%plot the early image
disp('Plotting the early image-DEM pair for the fjord');
figure1 = figure;
imagesc(A.x,A.y,A.z); hold on; 
colormap gray; set(gca,'clim',[1.05*min(A.z(~isnan(A.z))) (0.95)*max(max(A.z))]);
set(gca,'ydir','normal'); hold on;
set(gcf,'position',[50 100 600 600]);
figure2 = figure;
set(gcf,'position',[650 100 600 600]);
z = DEM1.z_elpsd_adjust+SL(i).initial.coreg_z*ones(size(DEM1.z_elpsd_adjust));
idx = nearestneighbour([A.x(1) A.x(end)],DEM1.x); idy = nearestneighbour([A.y(1) A.y(end)],DEM1.y);
DEM_x = DEM1.x(idx(1):idx(2)); DEM_y = DEM1.y(idy(1):idy(2));
DEM_z = z(idy(1):idy(2),idx(1):idx(2));
imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
set(gca,'ydir','normal','clim',[0 50]); hold on;

%extract pixel areas
DEM_pixel_area = abs(DEM_x(1)-DEM_x(2)).*abs(DEM_y(1)-DEM_y(2)); %square meters
im_pixel_area = abs(A.x(1)-A.x(2)).*abs(A.y(1)-A.y(2)); %square meters

% ----------extract initial area & size info----------
disp('For each iceberg, extract size info');
if ~exist([dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_shapes/'],'dir')
    mkdir([dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);
end
cd([dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);

%create density profiles
clear FAC;
FAC(1) = firnair.median; FAC(2) = firnair.median-firnair.uncert; FAC(3) = firnair.median+firnair.uncert; %estimate firn air content
ft = fittype('917-(917-a)*exp(-x/b)');[f,~] = fit(density_depths',density_levels',ft,'StartPoint',[350 firnair.median]); ci = confint(f); %create density profile
density_z = [0:1:1000];
density_profile = rho_i-(rho_i-f.a)*exp(-density_z/f.b); mindensity_profile = rho_i-(rho_i-ci(1,1))*exp(-density_z/ci(1,2)); maxdensity_profile = rho_i-(rho_i-ci(2,1))*exp(-density_z/ci(2,2));
%calculate wet density profile by flipping the shape of the exponential curve & compressing the range from (830-0) to (1026-830)
wetdensity_profile = 830+((830-density_profile)./830).*(rho_sw-830); wetdensity_profile(ceil(density.eightthir)+1:end) = density_profile(ceil(density.eightthir)+1:end);
minwetdensity_profile = 830+((830-mindensity_profile)./830).*(rho_sw-830); minwetdensity_profile(ceil(density.eightthir)+1:end) = mindensity_profile(ceil(density.eightthir)+1:end);
maxwetdensity_profile = 830+((830-maxdensity_profile)./830).*(rho_sw-830); maxwetdensity_profile(ceil(density.eightthir)+1:end) = maxdensity_profile(ceil(density.eightthir)+1:end);

%loop
for i = 1:size(SL,2)
    disp(['Extracting iceberg size info for iceberg #',num2str(i),' of ',num2str(size(SL,2))]);
    
    %zoom in
    disp(['...zooming-in & plotting the DEM ROI in the image and DEM for iceberg #',num2str(i)]);
    figure(figure1);
%     imagesc(A.x,A.y,A.z); set(gca,'ydir','normal'); hold on;
    colormap gray; set(gca,'clim',[1.05*min(A.z(~isnan(A.z))) (0.95)*max(max(A.z))]);
    vxi = nearestneighbour(SL(i).initial.x,A.x); vyi = nearestneighbour(SL(i).initial.y,A.y);
    set(gca,'xlim',[min(A.x(vxi))-150 max(A.x(vxi))+150],'ylim',[min(A.y(vyi))-150 max(A.y(vyi))+150]);
    prompt = 'Widen the zoom window (y/n)?';
    zstr = input(prompt,'s');
    if strmatch(zstr,'y')==1
        set(gca,'xlim',[min(A.x(vxi))-250 max(A.x(vxi))+250],'ylim',[min(A.y(vyi))-250 max(A.y(vyi))+250]);
    end
    plot(A.x(vxi),A.y(vyi),'--r','linewidth',2);
    figure(figure2);
%     imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    if strmatch(zstr,'y')==1
        set(gca,'xlim',[min(SL(i).initial.x)-250 max(SL(i).initial.x)+250],'ylim',[min(SL(i).initial.y)-250 max(SL(i).initial.y)+250]);
    else
        set(gca,'xlim',[min(SL(i).initial.x)-150 max(SL(i).initial.x)+150],'ylim',[min(SL(i).initial.y)-150 max(SL(i).initial.y)+150]);
    end
    plot(SL(i).initial.x,SL(i).initial.y,'-k','linewidth',2); hold on; plot(SL(i).initial.x,SL(i).initial.y,'--w','linewidth',2); hold on;
    
    %define the iceberg area & save its vertices to a shapefile
    disp('Draw a polygon following the iceberg edges in the WV image (used to estimate surface area).');
    figure(figure1);
    [iceberg_IMmask,x,y] = roipoly;
    
    %apply the mask from the image to the DEM
    [image_xgrid,image_ygrid] = meshgrid(A.x,A.y);
    [DEM_xgrid,DEM_ygrid] = meshgrid(DEM_x,DEM_y);
    iceberg_DEMmask = interp2(image_xgrid,image_ygrid,double(iceberg_IMmask),DEM_xgrid,DEM_ygrid);
    iceberg_DEMmask = round(iceberg_DEMmask);
    DEM_z_masked = iceberg_DEMmask.*DEM_z;
    DEM_z_masked(iceberg_DEMmask == 0) = 0;
    disp('NOTE: DEM coordinates may be shifted relative to the image: compare the DEM ROI (black & white dash) & the iceberg polygon (pink & white dash)on the DEM');
    figure(figure2);
    plot(x,y,'-w','linewidth',2); hold on; plot(x,y,'--m','linewidth',2); hold on;
    SL(i).initial.x_perim = x; SL(i).initial.y_perim = y;
    
    %modify the mask using the DEM if necessary
    disp('Now modify the iceberg mask to extract elevations from the entire iceberg area');
    clear DEM_z_masked iceberg_DEMmask x y;
    figure(figure2);
    [iceberg_DEMmask,x,y] = roipoly;
    SL(i).initial.x_perim = []; SL(i).initial.y_perim = [];
    SL(i).initial.x_perim = x; SL(i).initial.y_perim = y;
    DEM_z_masked = iceberg_DEMmask.*DEM_z;
    DEM_z_masked(iceberg_DEMmask == 0) = 0;
%     clf(figure2); drawnow;
%     imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
%     set(gca,'ydir','normal','clim',[0 50]); hold on;
%     set(gca,'xlim',[min(SL(i).initial.x)-200 max(SL(i).initial.x)+200],'ylim',[min(SL(i).initial.y)-200 max(SL(i).initial.y)+200]);
    plot(SL(i).initial.x,SL(i).initial.y,'-k','linewidth',2); hold on; plot(SL(i).initial.x,SL(i).initial.y,'--w','linewidth',2); hold on;
    plot(SL(i).initial.x_perim,SL(i).initial.y_perim,'-w','linewidth',2); hold on; plot(SL(i).initial.x_perim,SL(i).initial.y_perim,'--m','linewidth',2); hold on;
    
    %flag the iceberg as upright or overturned
    answer = questdlg('Is the iceberg upright (i.e., does it look like the glacier surface)?',...
    'Iceberg Upright or Flipped','1) Yes','2) No','1) Yes');
    switch answer
        case '1) Yes'
            SL(i).orientation = 1; %upright
        case '2) No'
            SL(i).orientation = 0; %overturned or a fragment of the full thickness
    end
    
    %use the mask to extract size & shape info
    disp('Pulling plan-view info');
    for k = 1:size(x,1)-1
        dist(k) = ((x(k)-x(k+1))^2+(y(k)-y(k+1))^2)^(1/2);
    end
    dist(end+1) = ((x(end)-x(1))^2+(y(end)-y(1))^2)^(1/2);
    shapes_o = regionprops(iceberg_IMmask,'Area','MinorAxisLength','MajorAxisLength'); %use the image-derived mask to estimate area
    lmin_o = abs(A.x(2)-A.x(1))*shapes_o.MinorAxisLength;
    lmax_o = abs(A.x(2)-A.x(1))*shapes_o.MajorAxisLength;
    SL(i).initial.width_min = lmin_o; SL(i).initial.width_max = lmax_o;
    base_area = shapes_o.Area.*im_pixel_area; SL(i).initial.SA = base_area;
    SL(i).initial.z_map = single(DEM_z_masked); SL(i).initial.z_map(SL(i).initial.z_map == 0) = NaN; SL(i).initial.z_map(SL(i).initial.z_map<0) = NaN;
    zo = SL(i).initial.z_map(~isnan(SL(i).initial.z_map));
    SL(i).initial.z_median = nanmedian(zo); SL(i).initial.z_mad = mad(zo,1);
    perimeter = sum(dist);
    SL(i).initial.radius = perimeter/(2*pi);
    clear iceberg_IMmask;
    
    %iteratively estimate bulk density
    disp('estimating density');
    if SL(i).orientation == 0 %berg is a fragment or overturned so assume the full firn is saturated
        for k = 1:3 %k=1 is the best guess, k=[2,3] constraints uncertainties
            j=1;
            Hberg(j) = (rho_sw/(rho_sw-rho_i))*SL(i).initial.z_median;
            Hberg_err(j) = abs(Hberg(j)).*sqrt((rho_sw_err/rho_sw)^2 + ((1.4826*SL(i).initial.z_mad)/SL(i).initial.z_median)^2 + ((rho_sw_err^2 + rho_i_err^2)/(rho_sw-rho_i)^2));
            wet_ref = find(density_z<=density.eightthir,1,'last');
            if Hberg(j) > density_z(wet_ref)
                rho_prof = [wetdensity_profile(1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [minwetdensity_profile(1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxwetdensity_profile(1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
            else
                rho_prof = [wetdensity_profile(1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [minwetdensity_profile(1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxwetdensity_profile(1:ceil(Hberg(j))+1)];
            end
            rho_f(j) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
            rho_f_err(j) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]),[],'omitnan');
            clear rho_prof*;
%             rho_f(j) = rho_i*((Hberg(j)-FAC(k)+((rho_sw/rho_i)*FAC(k)))./Hberg(j)); %assume firn is water-saturated
%             rho_f_err(j) = sqrt(rho_i_err^2 + ((rho_sw-rho_i)*FAC(k)/Hberg(j))^2*(((rho_sw_err^2+rho_i_err^2)/((rho_sw-rho_i)^2)) + ((1.4826*SL(i).initial.z_mad)/SL(i).initial.z_median)^2)); %excludes FAC uncertainty b/c I am solving for that directly
            while j
                Hberg(j+1) = (rho_sw/(rho_sw-rho_f(j)))*SL(i).initial.z_median;
                Hberg_err(j+1) = abs(Hberg(j+1)).*sqrt(((abs(rho_sw/(rho_sw-rho_f(j)))*sqrt((rho_sw_err/rho_sw)^2 + (sqrt(rho_sw_err^2+rho_f_err(j)^2)/(rho_sw-rho_f(j)))^2))/(rho_sw/(rho_sw-rho_f(j))))^2 + ((1.4826*SL(i).initial.z_mad)/SL(i).initial.z_median)^2);
                if Hberg(j+1) > density_z(wet_ref)
                    rho_prof = [wetdensity_profile(1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [minwetdensity_profile(1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxwetdensity_profile(1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                else
                    rho_prof = [wetdensity_profile(1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [minwetdensity_profile(1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxwetdensity_profile(1:ceil(Hberg(j+1))+1)];
                end
                rho_f(j+1) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
                rho_f_err(j+1) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]),[],'omitnan');
                clear rho_prof*;
                
                if abs(rho_f(j+1)-rho_f(j)) < 0.25*rho_f_err(j+1)
                    if k == 1
                        SL(i).density_type = 'wet';
                        SL(i).initial.density = rho_f(j+1);
                        %propagated density uncertainty
                        SL(i).initial.density_uncert = rho_f_err(j);
                    else
                        %density uncertainty using FAC uncertainty bounds to iterate thickness & rho_f estimates
                        SL(i).initial_range.density(k-1) = rho_f(j+1);
                    end
                    clear Hberg rho_f* dry_ref wet_ref;
                    break
                else
                    j = j+1;
                end
%                 clear Hice Vice;
            end
        end
    else %berg is upright so only wet the firn below the waterline
        for k = 1:3 %k=1 is the best guess, k=[2,3] constraints uncertainties
            j=1;
            Hberg(j) = (rho_sw/(rho_sw-rho_i))*SL(i).initial.z_median;
            Hberg_err(j) = abs(Hberg(j)).*sqrt((rho_sw_err/rho_sw)^2 + ((1.4826*SL(i).initial.z_mad)/SL(i).initial.z_median)^2 + ((rho_sw_err^2 + rho_i_err^2)/(rho_sw-rho_i)^2));
            dry_ref = find(density_z<=SL(i).initial.z_median,1,'last'); wet_ref = find(density_z<=density.eightthir,1,'last');
            if Hberg(j) > density_z(wet_ref)
                rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
            else
                rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:ceil(Hberg(j))+1)];
            end
            rho_f(j) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
            rho_f_err(j) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]),[],'omitnan');
            clear rho_prof*;

            while j
                Hberg(j+1) = (rho_sw/(rho_sw-rho_f(j)))*SL(i).initial.z_median;
                Hberg_err(j+1) = abs(Hberg(j+1)).*sqrt(((abs(rho_sw/(rho_sw-rho_f(j)))*sqrt((rho_sw_err/rho_sw)^2 + (sqrt(rho_sw_err^2+rho_f_err(j)^2)/(rho_sw-rho_f(j)))^2))/(rho_sw/(rho_sw-rho_f(j))))^2 + ((1.4826*SL(i).initial.z_mad)/SL(i).initial.z_median)^2);
                if Hberg(j+1) > density_z(wet_ref)
                    rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                else
                    rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:ceil(Hberg(j+1))+1)];
                end
                rho_f(j+1) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
                rho_f_err(j+1) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]),[],'omitnan');
                clear rho_prof*;
                
                if abs(rho_f(j+1)-rho_f(j)) < 0.25*rho_f_err(j+1)
                    if k == 1
                        SL(i).density_type = 'dry';
                        SL(i).initial.density = rho_f(j+1);
                        %propagated density uncertainty
                        SL(i).initial.density_uncert = rho_f_err(j);
                    else
                        %density uncertainty using FAC uncertainty bounds to iterate thickness & rho_f estimates
                        SL(i).initial_range.density(k-1) = rho_f(j+1);
                    end
                    clear Hberg rho_f* dry_ref wet_ref;
                    break
                else
                    j = j+1;
                end
%                 clear Hice Vice;
            end
        end
    end
    
    %estimate the iceberg depth & submerged area using densities from the Ligtenberg FDM
    disp('estimating thickness & volume');
    if SL(i).orientation == 0 %overturned/fragment
        rho_f = sort([rho_i SL(i).initial.density]); %figure out if bubble-free ice or water-saturated firn has a lower density
        for l = 1:2
            draft(l) = (rho_f(l)/(rho_sw-rho_f(l)))*SL(i).initial.z_median;
            Hberg(l) = (rho_sw/(rho_sw-rho_f(l)))*SL(i).initial.z_median;
        end
    else
        rho_f = sort(SL(i).initial_range.density); %assume density constrained by model FAC range
        for l = 1:2
            draft(l) = (rho_f(l)/(rho_sw-rho_f(l)))*SL(i).initial.z_median;
            Hberg(l) = (rho_sw/(rho_sw-rho_f(l)))*SL(i).initial.z_median;
        end
    end
    SL(i).initial_range.draft = sort(draft);
    lat_area = sort(draft).*sum(dist); SL(i).initial_range.LA = lat_area;
    area = base_area*ones(size(lat_area)) + lat_area; SL(i).initial_range.TA = area;
    SL(i).initial.V = DEM_pixel_area*(rho_sw/(rho_sw-SL(i).initial.density))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(i).initial.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(i).initial_range.V(1) = DEM_pixel_area*(rho_sw/(rho_sw-min(rho_f)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(i).initial.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(i).initial_range.V(2) = DEM_pixel_area*(rho_sw/(rho_sw-max(rho_f)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(i).initial.z_median*sum(sum(isnan(DEM_z_masked))));
    
    %save the iceberg outline
    disp('Writing shapefile containing early iceberg polygon information');
    S.Geometry = 'Polygon';
    S.BoundingBox = [min(x) min(y); max(x) max(y)];
    S.X = double(x'); S.Y = double(y');
    S.Name = ['iceberg',num2str(berg_numbers(i).name(8:9))];
    shapefile_name = ['WV_',num2str(to),'_icebergshape',num2str(berg_numbers(i).name(8:9))];
    shapewrite(S,shapefile_name);
    cd_to_site_data = ['cd ',dir_output]; eval(cd_to_site_data);
    copyfile([dir_code,'general/antarctic_PSprojection.prj'],[dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_shapes/',shapefile_name,'.prj']);
    cd_output_dir = ['cd ',dir_output,'/',DEM1.time,'-',DEM2.time,'/']; eval(cd_output_dir);
    clear S;

    %save the data
    disp('Saving data');
    save([dir_output,DEM1.time,'-',DEM2.time,'/',region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat'],'SL','-v7.3');
    cd([dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);
    
    if i<size(SL,2)
        disp(' and advancing to the next iceberg');
    end
%     clf(figure1); clf(figure2); drawnow;
    clear DEM_z_masked z idx idy x y iceberg_IMmask iceberg_DEMmask rho_f draft Hberg lat_area base_area dist area_mask masked_pixels image_xgrid image_ygrid shapes lmin* lmax* perimeter;
end
clear Z x1 x2 y1 y2 xlims ylims xmin xmax ymin ymax DEM_x DEM_y DEM_z ft f ci density_z *density_profile;
disp('Advancing to the later date');
close(figure1); close(figure2); clear A;

% ----------load the final image & DEM files to map the final iceberg margins----------

%set image bounds & crop
xf = []; yf = [];
for i = 1:length(SL)
    xf = [xf SL(i).final.x]; yf = [yf SL(i).final.y];
end
dy = IM2.y(1)-IM2.y(2);
if dy < 0
    x1 = find(IM2.x >= min(xf)-1500,1,'first'); x2 = find(IM2.x <= max(xf)+1500,1,'last');
    y1 = find(IM2.y >= max(yf)+1500,1,'first'); y2 = find(IM2.y <= min(yf)-1500,1,'last');
    xlims = sort([x1, x2]); ylims = sort([y1, y2]);
    xmin=xlims(1); xmax=xlims(2);
    xmin(isempty(xmin)) = 1; xmax(isempty(xmax)) = length(IM2.x);
    if length(ylims) ==1
        y2(isempty(y2)) = 1; y1(isempty(y1)) = length(IM2.y);
    end
    ylims = sort([y1, y2]);
    ymin=ylims(1); ymax=ylims(2);
else
    x1 = find(IM2.x >= min(xf)-1500,1,'first'); x2 = find(IM2.x <= max(xf)+1500,1,'last');
    y1 = find(IM2.y >= max(yf)+1500,1,'last'); y2 = find(IM2.y <= min(yf)-1500,1,'first');
    xlims = sort([x1, x2]); ylims = sort([y1, y2]);
    xmin=xlims(1); xmax=xlims(2);
    xmin(isempty(xmin)) = 1; xmax(isempty(xmax)) = length(IM2.x);
    if length(ylims)==1 
        y1(isempty(y1)) = 1; y2(isempty(y2)) = length(IM2.y);
    end
    ylims = sort([y1, y2]);
    ymin=ylims(1); ymax=ylims(2);
end
A.x = IM2.x(1,[xlims(1):xlims(2)]); A.y = IM2.y(1,[ylims(1):ylims(2)]);
A.z = double(IM2.z(ylims(1):ylims(2),xlims(1):xlims(2)));
if mean(A.z(~isnan(A.z)),'omitnan') < (1/3)*255 %image is dark
    z_adjust = imadjust(A.z./255,[],[],0.25);
    A.z = z_adjust; clear z_adjust;
end

%plot the early image
disp('Plotting the later image-DEM pair for the fjord');
figure1 = figure;
imagesc(A.x,A.y,A.z); hold on; 
colormap gray; set(gca,'clim',[1.05*min(A.z(~isnan(A.z))) (0.95)*max(max(A.z))]);
set(gca,'ydir','normal'); hold on;
set(gcf,'position',[50 100 600 600]);
figure2 = figure;
set(gcf,'position',[650 100 600 600]);
z = DEM2.z_elpsd_adjust+SL(i).final.coreg_z*ones(size(DEM2.z_elpsd_adjust));
idx = nearestneighbour([A.x(1) A.x(end)],DEM2.x); idy = nearestneighbour([A.y(1) A.y(end)],DEM2.y);
DEM_x = DEM2.x(idx(1):idx(2)); DEM_y = DEM2.y(idy(1):idy(2));
DEM_z = z(idy(1):idy(2),idx(1):idx(2));
imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
set(gca,'ydir','normal','clim',[0 50]); hold on;

%extract pixel areas
DEM_pixel_area = abs(DEM_x(1)-DEM_x(2)).*abs(DEM_y(1)-DEM_y(2)); %square meters
im_pixel_area = abs(A.x(1)-A.x(2)).*abs(A.y(1)-A.y(2)); %square meters

% ----------extract final area & size info----------
disp('For each iceberg, measure the iceberg aerial extent by drawing a polygon around the iceberg edge');
cd([dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);

%create density profiles
clear FAC;
FAC(1) = firnair.median; FAC(2) = firnair.median-firnair.uncert; FAC(3) = firnair.median+firnair.uncert; %estimate firn air content
ft = fittype('917-(917-a)*exp(-x/b)');[f,~] = fit(density_depths',density_levels',ft,'StartPoint',[350 firnair.median]); ci = confint(f); %create density profile
density_z = [0:1:1000];
density_profile = rho_i-(rho_i-f.a)*exp(-density_z/f.b); mindensity_profile = rho_i-(rho_i-ci(1,1))*exp(-density_z/ci(1,2)); maxdensity_profile = rho_i-(rho_i-ci(2,1))*exp(-density_z/ci(2,2));
%calculate wet density profile by flipping the shape of the exponential curve & compressing the range from (830-0) to (1026-830)
wetdensity_profile = 830+((830-density_profile)./830).*(rho_sw-830); wetdensity_profile(ceil(density.eightthir)+1:end) = density_profile(ceil(density.eightthir)+1:end);
minwetdensity_profile = 830+((830-mindensity_profile)./830).*(rho_sw-830); minwetdensity_profile(ceil(density.eightthir)+1:end) = mindensity_profile(ceil(density.eightthir)+1:end);
maxwetdensity_profile = 830+((830-maxdensity_profile)./830).*(rho_sw-830); maxwetdensity_profile(ceil(density.eightthir)+1:end) = maxdensity_profile(ceil(density.eightthir)+1:end);

%loop
for i = 1:size(SL,2) 
    disp(['Extracting iceberg size info for iceberg #',num2str(i),' of ',num2str(size(SL,2))]);
    
    %zoom in
    disp(['...zooming-in & plotting the DEM ROI in the image and DEM for iceberg #',num2str(i)]);
    figure(figure1);
    colormap gray; set(gca,'clim',[1.05*min(A.z(~isnan(A.z))) (0.95)*max(max(A.z))]);
    vxf = nearestneighbour(SL(i).final.x,A.x); vyf = nearestneighbour(SL(i).final.y,A.y);
    set(gca,'xlim',[min(A.x(vxf))-150 max(A.x(vxf))+150],'ylim',[min(A.y(vyf))-150 max(A.y(vyf))+150]);
    prompt = 'Widen the zoom window (y/n)?';
    zstr = input(prompt,'s');
    if strmatch(zstr,'y')==1
        set(gca,'xlim',[min(A.x(vxf))-250 max(A.x(vxf))+250],'ylim',[min(A.y(vyf))-250 max(A.y(vyf))+250]);
    end
    plot(A.x(vxf),A.y(vyf),'--r','linewidth',2);
    figure(figure2);
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    if strmatch(zstr,'y')==1
        set(gca,'xlim',[min(SL(i).final.x)-250 max(SL(i).final.x)+250],'ylim',[min(SL(i).final.y)-250 max(SL(i).final.y)+250]);
    else
        set(gca,'xlim',[min(SL(i).final.x)-150 max(SL(i).final.x)+150],'ylim',[min(SL(i).final.y)-150 max(SL(i).final.y)+150]);
    end
    plot(SL(i).final.x,SL(i).final.y,'-k','linewidth',2); hold on; plot(SL(i).final.x,SL(i).final.y,'--w','linewidth',2); hold on;
    
    %define the iceberg area & save its vertices to a shapefile
    disp('Draw a polygon following the iceberg edges in the WV image (used to estimate surface area).');
    figure(figure1);
    [iceberg_IMmask,x,y] = roipoly;
    
    %apply the mask from the image to the DEM
    [image_xgrid,image_ygrid] = meshgrid(A.x,A.y);
    [DEM_xgrid,DEM_ygrid] = meshgrid(DEM_x,DEM_y);
    iceberg_DEMmask = interp2(image_xgrid,image_ygrid,double(iceberg_IMmask),DEM_xgrid,DEM_ygrid);
    iceberg_DEMmask = round(iceberg_DEMmask);
    DEM_z_masked = iceberg_DEMmask.*DEM_z;
    DEM_z_masked(iceberg_DEMmask == 0) = 0;
    disp('NOTE: DEM coordinates may be shifted relative to the image: compare the DEM ROI (black & white dash) & the iceberg polygon (pink & white dash)on the DEM');
    figure(figure2);
    plot(x,y,'-w','linewidth',2); hold on; plot(x,y,'--m','linewidth',2); hold on;
    SL(i).final.x_perim = x; SL(i).final.y_perim = y;
    
    %modify the mask using the DEM if necessary
    disp('Now modify the iceberg mask need to extract elevations from the entire iceberg area');
    clear DEM_z_masked iceberg_DEMmask x y;
    figure(figure2);
    [iceberg_DEMmask,x,y] = roipoly;
    SL(i).final.x_perim = []; SL(i).final.y_perim = [];
    SL(i).final.x_perim = x; SL(i).final.y_perim = y;
    DEM_z_masked = iceberg_DEMmask.*DEM_z;
    DEM_z_masked(iceberg_DEMmask == 0) = 0;
    plot(SL(i).final.x,SL(i).final.y,'-k','linewidth',2); hold on; plot(SL(i).final.x,SL(i).final.y,'--w','linewidth',2); hold on;
    plot(SL(i).final.x_perim,SL(i).final.y_perim,'-w','linewidth',2); hold on; plot(SL(i).final.x_perim,SL(i).final.y_perim,'--m','linewidth',2); hold on;
    
    %use the mask to extract size & shape info
    disp('Pulling plan-view info');
    for k = 1:size(x,1)-1
        dist(k) = ((x(k)-x(k+1))^2+(y(k)-y(k+1))^2)^(1/2);
    end
    dist(end+1) = ((x(end)-x(1))^2+(y(end)-y(1))^2)^(1/2);
    shapes_f = regionprops(iceberg_IMmask,'Area','MinorAxisLength','MajorAxisLength');
    lmin_f = abs(A.x(2)-A.x(1))*shapes_f.MinorAxisLength;
    lmax_f = abs(A.x(2)-A.x(1))*shapes_f.MajorAxisLength;
    SL(i).final.width_min = lmin_f; SL(i).final.width_max = lmax_f;
    base_area = shapes_f.Area.*im_pixel_area; SL(i).final.SA = base_area;
    SL(i).final.z_map = single(DEM_z_masked); SL(i).final.z_map(SL(i).final.z_map == 0) = NaN; SL(i).final.z_map(SL(i).final.z_map<0) = NaN;
    zf = SL(i).final.z_map(~isnan(SL(i).final.z_map));
    SL(i).final.z_median = nanmedian(zf); SL(i).final.z_mad = mad(zf,1);
    perimeter = sum(dist);
    SL(i).final.radius = perimeter/(2*pi);
    clear iceberg_IMmask;
    
    %iteratively estimate bulk density
    disp('estimating density');
    if SL(i).orientation == 0 %berg is a fragment or overturned so assume firn is totally saturated
        for k = 1:3 %k=1 is the best guess, k=[2,3] constraints uncertainties
            j=1;
            Hberg(j) = (rho_sw/(rho_sw-rho_i))*SL(i).final.z_median;
            Hberg_err(j) = abs(Hberg(j)).*sqrt((rho_sw_err/rho_sw)^2 + ((1.4826*SL(i).final.z_mad)/SL(i).final.z_median)^2 + ((rho_sw_err^2 + rho_i_err^2)/(rho_sw-rho_i)^2));
            wet_ref = find(density_z<=density.eightthir,1,'last');
            if Hberg(j) > density_z(wet_ref)
                rho_prof = [wetdensity_profile(1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [minwetdensity_profile(1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxwetdensity_profile(1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
            else
                rho_prof = [wetdensity_profile(1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [minwetdensity_profile(1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxwetdensity_profile(1:ceil(Hberg(j))+1)];
            end
            rho_f(j) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
            rho_f_err(j) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]));
            clear rho_prof*;
            
            while j
                Hberg(j+1) = (rho_sw/(rho_sw-rho_f(j)))*SL(i).final.z_median;
                Hberg_err(j+1) = abs(Hberg(j+1)).*sqrt(((abs(rho_sw/(rho_sw-rho_f(j)))*sqrt((rho_sw_err/rho_sw)^2 + (sqrt(rho_sw_err^2+rho_f_err(j)^2)/(rho_sw-rho_f(j)))^2))/(rho_sw/(rho_sw-rho_f(j))))^2 + ((1.4826*SL(i).final.z_mad)/SL(i).final.z_median)^2);
                if Hberg(j+1) > density_z(wet_ref)
                    rho_prof = [wetdensity_profile(1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [minwetdensity_profile(1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxwetdensity_profile(1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                else
                    rho_prof = [wetdensity_profile(1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [minwetdensity_profile(1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxwetdensity_profile(1:ceil(Hberg(j+1))+1)];
                end
                rho_f(j+1) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
                rho_f_err(j+1) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]));
                clear rho_prof*;
                
                if abs(rho_f(j+1)-rho_f(j)) < 0.25*rho_f_err(j+1)
                    if k == 1
                        SL(i).density_type = 'wet';
                        SL(i).final.density = rho_f(j+1);
                        %propagated density uncertainty
                        SL(i).final.density_uncert = rho_f_err(j);
                    else
                        %density uncertainty using FAC uncertainty bounds to iterate thickness & rho_f estimates
                        SL(i).final_range.density(k-1) = rho_f(j+1);
                    end
                    clear Hberg rho_f* dry_ref wet_ref;
                    break
                else
                    j = j+1;
                end
%                 clear Hice Vice;
            end
        end
    else %berg is upright so only wet the firn below the waterline
        for k = 1:3 %k=1 is the best guess, k=[2,3] constraints uncertainties
            j=1;
            Hberg(j) = (rho_sw/(rho_sw-rho_i))*SL(i).final.z_median;
            Hberg_err(j) = abs(Hberg(j)).*sqrt((rho_sw_err/rho_sw)^2 + ((1.4826*SL(i).final.z_mad)/SL(i).final.z_median)^2 + ((rho_sw_err^2 + rho_i_err^2)/(rho_sw-rho_i)^2));
            dry_ref = find(density_z<=SL(i).final.z_median,1,'last'); wet_ref = find(density_z<=density.eightthir,1,'last');
            if Hberg(j) > density_z(wet_ref)
                rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
            else
                rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:ceil(Hberg(j))+1)];
            end
            rho_f(j) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
            rho_f_err(j) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]));
            clear rho_prof*;
            
            while j
                Hberg(j+1) = (rho_sw/(rho_sw-rho_f(j)))*SL(i).final.z_median;
                Hberg_err(j+1) = abs(Hberg(j+1)).*sqrt(((abs(rho_sw/(rho_sw-rho_f(j)))*sqrt((rho_sw_err/rho_sw)^2 + (sqrt(rho_sw_err^2+rho_f_err(j)^2)/(rho_sw-rho_f(j)))^2))/(rho_sw/(rho_sw-rho_f(j))))^2 + ((1.4826*SL(i).final.z_mad)/SL(i).final.z_median)^2);
                if Hberg(j+1) > density_z(wet_ref)
                    rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                else
                    rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:ceil(Hberg(j+1))+1)];
                end
                rho_f(j+1) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
                rho_f_err(j+1) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]));
                clear rho_prof*;
                
                if abs(rho_f(j+1)-rho_f(j)) < 0.25*rho_f_err(j+1)
                    if k == 1
                        SL(i).density_type = 'dry';
                        SL(i).final.density = rho_f(j+1);
                        %propagated density uncertainty
                        SL(i).final.density_uncert = rho_f_err(j);
                    else
                        %density uncertainty using FAC uncertainty bounds to iterate thickness & rho_f estimates
                        SL(i).final_range.density(k-1) = rho_f(j+1);
                    end
                    clear Hberg rho_f* dry_ref wet_ref;
                    break
                else
                    j = j+1;
                end
%                 clear Hice Vice;
            end
        end
    end
    
    %estimate the iceberg depth & submerged area using densities from the Ligtenberg FDM
    disp('Estimating thickness & volume');
    if SL(i).orientation == 0 %overturned/fragment
        rho_f = sort([rho_i SL(i).final.density]); %figure out if bubble-free ice or water-saturated firn has a lower density
        for l = 1:2
            draft(l) = (rho_f(l)/(rho_sw-rho_f(l)))*SL(i).final.z_median;
            Hberg(l) = (rho_sw/(rho_sw-rho_f(l)))*SL(i).final.z_median;
        end
    else
        rho_f = sort(SL(i).final_range.density); %assume density constrained by model FAC range
        for l = 1:2
            draft(l) = (rho_f(l)/(rho_sw-rho_f(l)))*SL(i).final.z_median;
            Hberg(l) = (rho_sw/(rho_sw-rho_f(l)))*SL(i).final.z_median;
        end
    end
    SL(i).final_range.draft = sort(draft);
    lat_area = sort(draft).*sum(dist); SL(i).final_range.LA = lat_area;
    area = base_area*ones(size(lat_area)) + lat_area; SL(i).final_range.TA = area;
    SL(i).final.V = DEM_pixel_area*(rho_sw/(rho_sw-SL(i).final.density))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(i).final.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(i).final_range.V(1) = DEM_pixel_area*(rho_sw/(rho_sw-min(rho_f)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(i).final.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(i).final_range.V(2) = DEM_pixel_area*(rho_sw/(rho_sw-max(rho_f)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(i).final.z_median*sum(sum(isnan(DEM_z_masked))));
    
    %save the iceberg outline
    S.Geometry = 'Polygon';
    S.BoundingBox = [min(x) min(y); max(x) max(y)];
    S.X = double(x'); S.Y = double(y');
    S.Name = ['iceberg',num2str(berg_numbers(i).name(8:9))];
    shapefile_name = ['WV_',num2str(tf),'_icebergshape',num2str(berg_numbers(i).name(8:9))];
    shapewrite(S,shapefile_name);
    copyfile([dir_code,'general/antarctic_PSprojection.prj'],[dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/',shapefile_name,'.prj']);
    cd([dir_output,DEM1.time,'-',DEM2.time,'/']);
    clear S;

    %calculate mean surface area for the exposed face (SA) &
    %submerged (TA) portions of the iceberg & associated temporal uncertainty
    SL(i).mean.z = (SL(i).initial.z_median+SL(i).final.z_median)/2;
    SL(i).mean.SA = (SL(i).initial.SA+SL(i).final.SA)/2;
    SL(i).change.SA = abs(SL(i).initial.SA-SL(i).final.SA)/2;
    SL(i).mean.TA = (nanmean(SL(i).initial_range.TA)+nanmean(SL(i).final_range.TA))/2;
    SL(i).uncert.TA = max([(SL(i).initial_range.TA(2)-SL(i).initial_range.TA(1))/2 (SL(i).final_range.TA(2)-SL(i).final_range.TA(1))/2]);
    SL(i).change.TA = abs(nanmean(SL(i).initial_range.TA)-nanmean(SL(i).final_range.TA))/2;
    SL(i).mean.V = (nanmean(SL(i).initial_range.V)+nanmean(SL(i).final_range.V))/2;
    SL(i).uncert.V = max([(SL(i).initial_range.V(2)-SL(i).initial_range.V(1))/2 (SL(i).final_range.V(2)-SL(i).final_range.V(1))/2]);
    SL(i).change.V = abs(nanmean(SL(i).initial_range.V)-nanmean(SL(i).final_range.V))/2;
    
    %save the data
    disp('Saving data');
    save([dir_output,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat'],'SL','-v7.3');
    cd([dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);
    
    if i<size(SL,2)
        disp(' and advancing to the next iceberg');
    end
%     clf(figure1); clf(figure2); drawnow;
    clear DEM_z_masked z idx idy x y iceberg_IMmask iceberg_DEMmask rho_f draft Hberg lat_area base_area dist area_mask masked_pixels image_xgrid image_ygrid shapes lmin* lmax* perimeter;
end
clear Y x1 x2 y1 y2 xlims ylims xmin xmax ymin ymax DEM_x DEM_y DEM_z ft f ci density_z *density_profile;
close(figure1); close(figure2); clear A;
   
% ----------use the iceberg surface area measurements to calculate volume fluxes & melt rates----------
disp('Converting changes in iceberg elevation to volume fluxes & melt rates');
close all; drawnow;
for i = 1:size(berg_numbers,1)
    %cd to the melt rate file & load each file sequentially
    cd(dir_iceberg);
    berg_number = berg_numbers(i).name(1:end-7);
    load(berg_numbers(i).name);
    
    %calculate the ice thickness (change to H = 3*(V/SA) if assuming the submerged shape is a cone!)
    SL(i).initial_range.H = (repmat(rho_sw,size(SL(i).initial_range.density))./(repmat(rho_sw,size(SL(i).initial_range.density))-SL(i).initial_range.density)).*SL(i).initial.z_median; 
    SL(i).final_range.H = (repmat(rho_sw,size(SL(i).final_range.density))./(repmat(rho_sw,size(SL(i).final_range.density))-SL(i).final_range.density)).*SL(i).final.z_median; 
    SL(i).mean.H = nanmean([SL(i).initial_range.H SL(i).final_range.H]); 
    SL(i).change.H = nanmean(SL(i).initial_range.H) - nanmean(SL(i).final_range.H);
    SL(i).mean.draft = nanmean([SL(i).initial_range.draft SL(i).final_range.draft]);
    SL(i).change.draft = nanmean(SL(i).initial_range.draft) - nanmean(SL(i).final_range.draft);
    
    %convert surface elevation changes to thickness changes
    vals = [];
    for j = 1:length(IB)
        if IB(j).local_adjust_f ~= 0 %use local sea level-adjusted elevations
            non_NaNs = length(IB(j).dz.local_adjust.map(~isnan(IB(j).dz.local_adjust.map)));
            vals(size(vals,2)+1:size(vals,2)+non_NaNs) = IB(j).dz.local_adjust.map(~isnan(IB(j).dz.local_adjust.map));
        else %use bedrock & tide-adjusted elevations
            non_NaNs = length(IB(j).dz.br_tide_adjust.map(~isnan(IB(j).dz.br_tide_adjust.map)));
            vals(size(vals,2)+1:size(vals,2)+non_NaNs) = IB(j).dz.br_tide_adjust.map(~isnan(IB(j).dz.br_tide_adjust.map));
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
    SL(i).mean.dz = nanmean(vals); SL(i).uncert.dz = nanstd(vals);
    dZ_mean = (rho_sw/(rho_sw-nanmean([SL(i).initial.density SL(i).final.density])))*SL(i).mean.dz; 
    dZ_stdev = (rho_sw/(rho_sw-nanmean([SL(i).initial.density SL(i).final.density])))*SL(i).uncert.dz;
    
    %correct for SMB loss
    dH_SMBadjust_mean = dZ_mean + SL(i).SMB; %SL(i).SMB=-surfmelt
    
    %correct for creep thinning
    rf_o = 3.5e-25; %rate factor at -10C (Pa^-3 s^-1) from Cuffey & Paterson p. 74
    if iceberg_avgtemp <= 263
        rf=rf_o*exp((-60000/8.314)*((1/SL(i).airtemp)-(1/263))); %rate factor for cold ice
    else
        rf=rf_o*exp((-134000/8.314)*((1/SL(i).airtemp)-(1/263))); %rate factor for nearly-temperate ice
    end
    SL(i).ratefactor = rf;
    B = rf^(-1/3); %Pa s^1/3
    creep = ((-1/(2*sqrt(3)))*((nanmean([SL(i).initial.density SL(i).final.density])*9.81*SL(i).mean.z)/(2*sqrt(3)))^3*(1-(nanmean([SL(i).initial.density SL(i).final.density])/rho_sw))^3)/(B^3); %creep thinning rate (1/s)
    SL(i).creep_dz = (SL(i).mean.H*creep*(dt*31536000));
    dH_submelt = dH_SMBadjust_mean + SL(i).creep_dz; %integrate creep over the ice thickness & over the time period

    %dH uncertainty sources
    %S1) Bulk density (use range in density estimates, including ice density if overturned)
    %S2) Ocean water density (use 1026 +/- 2 kg/m^3)
    %S3) Runoff uncertainty (from RACMO output)
    %R1) DEM uncertainty (random errors estimated as 2.9m)
    %R2) Translocation & Rotation (User) error (quantify by repeating procedure 10 times for each glacier)
    
    %%quantify potential bias from systematic errors
    %uncertainty from S1 & S2
    if SL(i).orientation == 0
        rho_f = nanmean([SL(i).initial.density SL(i).final.density]);
        rho_f_range = [rho_i SL(i).initial_range.density SL(i).final_range.density];
        rho_f_err = [min(rho_f_range-rho_f) max(rho_f_range-rho_f)];
    else
        rho_f = nanmean([SL(i).initial.density SL(i).final.density]);
        rho_f_range = [SL(i).initial_range.density SL(i).final_range.density];
        rho_f_err = [min(rho_f_range-rho_f) max(rho_f_range-rho_f)];
    end
    rho_added_err = sqrt(rho_f_err.^2 + repmat(rho_sw_err,size(rho_f_err)).^2);
    rho_conversion = rho_sw/(rho_sw-rho_f); 
    rho_conversion_err = abs(rho_conversion)*sqrt(repmat((rho_sw_err/rho_sw),size(rho_f_err)).^2 + (rho_added_err./repmat(rho_sw-rho_f,size(rho_f_err))).^2);
    
    %uncertainty from S3
    dz_SMB_err = 0.3*abs(surfmelt);
    
    %uncertainty from R1
    pixelno_o = SL(i).initial.SA./DEM_pixel_area; pixelno_f = SL(i).final.SA./DEM_pixel_area;
    dz_randerr = sqrt((2.9^2+2.9^2)./nanmean([pixelno_o pixelno_f])); %(2.9^2+2.9^2) is to account for DEM differencing over time

    %uncertainty from R2
    for k = 1:10
        dzs(k) = IB(k).dz.local_adjust.mean;
    end
    dz_user_err = nanstd(dzs);
    
    %total, systematic, and random errors
    SL(i).uncert.h_SMB = dz_SMB_err;
    SL(i).uncert.h_rand = dz_randerr;
    SL(i).uncert.h_user = dz_user_err;
    SL(i).uncert.hH_conversion = rho_conversion_err;
    dH_total_err = abs(dH_submelt).*sqrt(repmat(((dz_SMB_err^2 + dz_randerr^2 + dz_user_err^2)/SL(i).mean.dz^2),size(rho_conversion_err)) + (rho_conversion_err./repmat(rho_conversion,size(rho_conversion_err))).^2);
    
    %convert all 'dHs', which are change in thickness as estimated using the elevation change,
    %(& errors) to changes in volume so that the iceberg melt rate
    %(dH/dt) can be calculated by dividing the volume change by the total
    %submerged iceberg area
    dV_mean = SL(i).mean.SA*dH_submelt;
    dV_stdev = abs(dV_mean).*sqrt((dZ_stdev./dH_submelt)^2);
    dV_total_err = abs(dV_mean).*sqrt((dH_total_err./repmat(dH_submelt,size(dH_total_err))).^2 + repmat((SL(i).change.SA/2)./SL(i).mean.SA,size(dH_total_err)).^2);
    
    %convert to a melt rate (volume change with time & thickness change
    %over the entire submerged face)
    dVdt_mean = dV_mean/dt;  dVdt_stdev = dV_stdev/dt;
    dVdt_total_err = dV_total_err./dt;
    dHdt_mean = dVdt_mean/SL(i).mean.TA;
    dHdt_stdev = abs(dHdt_mean).*sqrt((dVdt_stdev./dVdt_mean)^2+(SL(i).uncert.TA./SL(i).mean.TA)^2);
    dHdt_total_err = abs(dHdt_mean).*sqrt((dVdt_total_err./repmat(dVdt_mean,size(dVdt_total_err))).^2 + repmat(SL(i).uncert.TA./SL(i).mean.TA,size(dVdt_total_err)).^2);
    
    %save volume fluxes and melt rates to the structure
    SL(i).mean.dVdt = dVdt_mean; 
    SL(i).range.dVdt(1) = dVdt_mean-3*dVdt_stdev; SL(i).range.dVdt(2) = dVdt_mean+3*dVdt_stdev; 
    SL(i).uncert.dVdt = dVdt_total_err;
    SL(i).mean.dHdt = dHdt_mean; 
    SL(i).range.dHdt(1) = dHdt_mean-3*dHdt_stdev; SL(i).range.dHdt(2) = dHdt_mean+3*dHdt_stdev; 
    SL(i).uncert.dHdt = dHdt_total_err;
    
%     %if overlapping bedrock regions were present, add the uniform offset that
%     %would be applied from bedrock & tidal differences to the structure
%     if ~isnan(IB(1).dz.br_tide_adjust.mean)
%         SL(i).elpsd_adjust_o = IB(1).elpsd_adjust_o; SL(i).elpsd_bt_adjust_f = IB(1).elpsd_adjust_f+IB(1).dbedrock_f+IB(1).dtide_f;
%     end
    
    
end

%save the compiled data
disp('Saving melt rates');
save([dir_iceberg,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat'],'SL','-v7.3');
disp('Melt rates saved! Advance to the next dataset');

close all; drawnow;

end

