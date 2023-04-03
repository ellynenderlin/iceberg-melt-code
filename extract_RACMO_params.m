function [dt,berg_T,berg_runoff,firnair,firn_density,firn_fit,firn_fitci] = extract_RACMO_params(dir_SMB,geography,berg_x,berg_y,berg_dates)

%read the RACMO files & find nearest pixel
cd(dir_SMB);
if geography == 0
    error('Use MAR data for Greenland!');
elseif geography == 1
    %check if the data fall within the RACMO AP model domain (preferred
    %over the AIS model b/c the spatial resolutions are 5.5 km & 27 km respectively)
    [berg_lon,berg_lat] = ps2wgs(berg_x,berg_y,'StandardParallel',-71,'StandardMeridian',0);
    runoff_lat = ncread('RACMO2.3p2_XPEN055_runoff_daily_2011_2016.nc','lat');
    runoff_lon = ncread('RACMO2.3p2_XPEN055_runoff_daily_2011_2016.nc','lon');
    for i = 1:size(runoff_lat,1)
        for j = 1:size(runoff_lat,2)
            [runoff_x(i,j),runoff_y(i,j)] = wgs2ps(runoff_lon(i,j),runoff_lat(i,j),'StandardParallel',-71,'StandardMeridian',0);
        end
    end
    in = inpolygon(berg_x,berg_y,[runoff_x(1,1) runoff_x(1,end) runoff_x(end,end) runoff_x(end,1) runoff_x(1,1)],...
        [runoff_y(1,1) runoff_y(1,end) runoff_y(end,end) runoff_y(end,1) runoff_y(1,1)]);
    %if the site is in RACMO AP, use those data
    if sum(in) >= 1
        AP=1;
        %load SMB params
        runoff = ncread('RACMO2.3p2_XPEN055_runoff_daily_2011_2016.nc','runoff'); runoff(runoff==-9999)=NaN;
        snowmelt = ncread('RACMO2.3p2_XPEN055_snowmelt_daily_2011_2016.nc','snowmelt'); snowmelt(snowmelt==-9999)=NaN;
        airtemp = ncread('RACMO2.3p2_XPEN055_T2m_daily_2011_2016.nc','t2m'); airtemp(airtemp==-9999)=NaN;
        icetemp = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','t2m'); icetemp(icetemp==-9999)=NaN; %assume the ice temp matches the avg long-term temp (for creep estimates)
        smb = [];
        runoff_time = ncread('RACMO2.3p2_XPEN055_runoff_daily_2011_2016.nc','time');
    else %otherwise use RACMO Antarctica
        clear runoff_lat runoff_lon runoff_x runoff_y;
        AP=0;
        runoff_lat = ncread('RACMO2.3p2_ANT27_runoff_daily_2011_2016.nc','lat');
        runoff_lon = ncread('RACMO2.3p2_ANT27_runoff_daily_2011_2016.nc','lon');
        for i = 1:size(runoff_lat,1)
            for j = 1:size(runoff_lat,2)
                [runoff_x(i,j),runoff_y(i,j)] = wgs2ps(runoff_lon(i,j),runoff_lat(i,j),'StandardParallel',-71,'StandardMeridian',0);
            end
        end
        %load SMB params
        runoff = ncread('RACMO2.3p2_ANT27_runoff_daily_2011_2016.nc','runoff'); runoff(runoff==-9999)=NaN;
        snowmelt = [];
        airtemp = ncread('RACMO2.3p2_ANT27_T2m_daily_2011_2016.nc','t2m'); airtemp(airtemp==-9999)=NaN;
        icetemp = ncread('RACMO2.3p2_ANT27_T2m_monthly_1979_2016.nc','t2m'); icetemp(icetemp==-9999)=NaN; %assume the ice temp matches the avg long-term temp (for creep estimates)
        smb = ncread('RACMO2.3p2_ANT27_smb_daily_2011_2016.nc','smb'); smb(smb==-9999)=NaN;
        runoff_time = ncread('RACMO2.3p2_ANT27_runoff_daily_2011_2016.nc','time');
        
    end
    clear in;
end
runoff_days = runoff_time-runoff_time(1); %days since 20110101

%calculate the x&y distance between the target of interest (berg lon, berg lat) and each RACMO grid cell
lat_diff = abs(berg_lat - runoff_lat);
lon_diff = abs(berg_lon - runoff_lon);
diff_map = sqrt(lat_diff.^2+lon_diff.^2); diff_map(isnan(squeeze(nanmean(runoff(:,:,1,270:450),4)))) = NaN;%solve for the distance vector using the x&y distances
RACMO_ref = find(diff_map==min(min(diff_map))); %find the minimum distance (reference for your grid cell is output)
[RACMOy RACMOx] = ind2sub(size(squeeze(nanmean(runoff(:,:,1,270:450),4))),RACMO_ref); %convert cell reference to an x- and y-cell index
disp(['RACMO x-reference = ',num2str(RACMOx),' & y-reference = ',num2str(RACMOy)]);

%adjust RACMO reference grid cell if necessary
disp('Adjust coordinates (if necessary) to extract surface melt estimates');
figure; set(gcf,'position',[100 100 700 700]);
runoff_cmap = colormap(jet(10001)); runoff_cmap(1,:) = [1 1 1];
if ~isempty(snowmelt)
    snowmelt_map = max(snowmelt(:,:,270:450),[],3).*86500./1000; snowmelt_map(isnan(snowmelt_map)) = min(snowmelt_map(~isnan(snowmelt_map)))-1;
    imagesc(snowmelt_map); colormap(gca,runoff_cmap); hold on; set(gca,'ydir','reverse'); 
    disp(['snowmelt at RACMO pixel = ',num2str(max(snowmelt(RACMOy,RACMOx,270:450)).*86500./1000),' mm w.e.']);
else
    smb_map = max(smb(:,:,270:450),[],3).*86500./1000; smb_map(isnan(smb_map)) = min(smb_map(~isnan(smb_map)))-1;
    imagesc(smb_map); colormap(gca,runoff_cmap); hold on; set(gca,'ydir','reverse'); 
    disp(['surface mass balance at RACMO pixel = ',num2str(max(smb(RACMOy,RACMOx,270:450)).*86500./1000),' mm w.e.']);
end
set(gca,'clim',[0 0.05]); cbar = colorbar; set(get(cbar,'ylabel'),'string','Runoff (m w.e. per day');
plot(RACMOx,RACMOy,'ok','markerfacecolor','k','markeredgecolor','w','markersize',20); hold on;
plot(RACMOx,RACMOy,'xw','markersize',20); hold on;
set(gca,'xlim',[RACMOx-10 RACMOx+10],'ylim',[RACMOy-10 RACMOy+10]);

%adjust coordinates as necessary
disp('Modify coordinates if the marker is in a region with no data');
adjustment = questdlg('Do the coordinates need to be modified?',...
    'RACMO coordinate adjustment','yes','no','yes');
switch adjustment
    case 'yes'
        disp('To change coordinates, identify new coordinates using tick marks and type RACMOx=XX; RACMOy=YY; dbcont ');
        keyboard
    case 'no'
        disp('Using default coordinates');
end
clear adjustment;
close all; drawnow;

%calculate the time separation between DEMs in terms of
%decimal years (ddays) & decimal days (days)
to = berg_dates(1,:); tf = berg_dates(2,:);
dt = datenum(tf(1:12),'yyyymmddHHMM') - datenum(to(1:12),'yyyymmddHHMM');


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
if mod(str2num(to(1:4)),4)==0; modayso = [31 29 31 30 31 30 31 31 30 31 30 31]; else modayso = [31 28 31 30 31 30 31 31 30 31 30 31]; end
if mod(str2num(tf(1:4)),4)==0; modaysf = [31 29 31 30 31 30 31 31 30 31 30 31]; else modaysf = [31 28 31 30 31 30 31 31 30 31 30 31]; end
doyo = sum(modayso(1:str2num(to(5:6))))-31+str2num(to(7:8)); doyf = sum(modaysf(1:str2num(tf(5:6))))-31+str2num(tf(7:8));
doy1 = sum(doys(1:str2num(to(1:4))-2011))+doyo;
doy2 = sum(doys(1:str2num(tf(1:4))-2011))+doyf;
decidayo =  str2num(to(1:4))+(doyo/doys(k)); decidayf =  str2num(tf(1:4))+(doyf/doys(k));
if decidayf<(2017+(max(runoff_days)-365-366-365-365-365-366)/365)
    melt = squeeze(runoff(RACMOy,RACMOx,doy1:doy2));
else
    if doy1 <= length(runoff)
        melt = cat(1,squeeze(runoff(RACMOy,RACMOx,doy1:length(runoff))),squeeze(runoff(RACMOy,RACMOx,sum(doys(1:2016-2011)):sum(doys(1:2016-2011))+doyf-1))); %temporarily use 2016 melt estimates from the same DOY until RACMO data are updated
    else
        if str2num(tf(1:4)) - str2num(to(1:4)) <= 1
            melt = cat(1,squeeze(runoff(RACMOy,RACMOx,sum(doys(1:2016-2011))+doyo:length(runoff))),squeeze(runoff(RACMOy,RACMOx,sum(doys(1:2016-2011)):sum(doys(1:2016-2011))+doyf-1)));
        else
            melt = cat(1,squeeze(runoff(RACMOy,RACMOx,sum(doys(1:2016-2011))+doyo:length(runoff))),squeeze(runoff(RACMOy,RACMOx,sum(doys(1:2016-2011)):sum(doys(1:2017-2011)))),squeeze(runoff(RACMOy,RACMOx,sum(doys(1:2016-2011)):sum(doys(1:2016-2011))+doyf-1)));
        end
    end
end
days = ones(1,length(melt)); %days(2:end-1) = 1; 
days(1) = ceil(datenum(to,'yyyymmddHHMMSS'))-datenum(to,'yyyymmddHHMMSS'); 
days(end) = datenum(tf,'yyyymmddHHMMSS')-floor(datenum(tf,'yyyymmddHHMMSS'));
% if length(days)~=length(melt)
%     clear days;
%     days = ones(1,length(melt)); days(1) = 1-hrs_o; days(2:end-1) = 1; days(end) = hrs_f;
% end
berg_runoff = nansum(days'.*melt)/1000; %surface meltwater that runs off (convert from mm w.e. per day to total m w.e.)
%estimate the ice temperature as the average long-term air temperature (doesn't account for advection of colder ice and melt/refreezing at the surface and/or submarine interface)
berg_T = nanmean(squeeze(icetemp(RACMOy,RACMOx,:))); % air temp (Kelvin)


%extract firn density estimates if in Antarctica, assume ice in Greenland
if geography == 1
    %extract FAC data from the Ligtenberg model output
    cd(dir_SMB); cd ../FDM_Antarctica/;
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
    figure; set(gcf,'position',[100 100 700 700]); 
    firn_cmap = colormap(jet(10001)); firn_cmap(1,:) = [1 1 1];
    firn_map = nanmean(FAC,3); firn_map(isnan(firn_map)) = min(firn_map(~isnan(firn_map))) - 1;
    imagesc(firn_map); colormap(gca,firn_cmap); hold on; set(gca,'ydir','reverse');
    cbar = colorbar; set(get(cbar,'ylabel'),'string','FAC (m) ');
    plot(firnx,firny,'ok','markerfacecolor','k','markeredgecolor','w','markersize',20); hold on;
    plot(firnx,firny,'xw','markersize',20); hold on;
    set(gca,'xlim',[firnx-10 firnx+10],'ylim',[firny-10 firny+10]);
    disp(['firn air content estimate = ',num2str(nanmean(FAC(firny,firnx,:),3)),'m']);
    
    %adjust coordinates as necessary
    disp('Modify coordinates if the marker is in a region with no data');
    adjustment = questdlg('Do the coordinates need to be modified?',...
        'FDM coordinate adjustment','yes','no','yes');
    switch adjustment
        case 'yes'
            disp('To change coordinates, identify new coordinates using tick marks and type firnx=XX; firny=YY; dbcont ');
            keyboard
        case 'no'
            disp('Using default coordinates');
    end
    clear adjustment;
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
    firn_density.xref = densityx; firn_density.yref = densityy;
    firn_density.sevhun = sevenhundred(firn_density.yref,firn_density.xref);
    firn_density.sevfif = sevenfifty(firn_density.yref,firn_density.xref);
    firn_density.eighthun = eighthundred(firn_density.yref,firn_density.xref);
    firn_density.eightthir = eightthirty(firn_density.yref,firn_density.xref);
    density_levels = [700 750 800 830]; density_depths = [firn_density.sevhun firn_density.sevfif firn_density.eighthun firn_density.eightthir];
    %fit a curve to the density-depth profile to estimate the depth of the base of the firn column
    [firn_fit,~] = firndensity_curvefit(density_depths,density_levels,firnair.median); firn_fitci = confint(firn_fit); %empirical, exponential density-depth relation from Schytt (1958) on Cuffey & Paterson p. 19
    
else %assume a uniform density for Greenland
    firnair = 0; firn_density = 900; firn_fit = NaN; firn_fitci = NaN;
end




end