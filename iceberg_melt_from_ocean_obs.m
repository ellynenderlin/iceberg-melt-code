%%% Extract ocean info from CTD and APB data around Antarctica

%%inputs:
%1) tab-delimited text files of iceberg melt data created by estimate_iceberg_melt_volume.m (*iceberg_meltinfo.txt)
%2) background image for Antarctic maps (Antarctic_RAMP_image_v2_1km.tif)
%3) ship-based conductivity, temperature, & depth mat-files (*CTD.mat) & autonomous pinneped data (*APB.mat): from Carlos Moffat
%4) RACMO data: RACMO2.3p2_XPEN055_T2m_daily_2011_2016.nc & RACMO2.3p2_ANT27_T2m_daily_2011_2016.nc

%%outputs:
%1) Antarctic-ocean-data.mat: CTD & APB structures contain compiled data in a standardized format
%2) Antarctic-icebergmelt-comparison.mat: melt structure w/ iceberg data & nearby CTD and/or APB data
%3) Antarctic-iceberg-oceandata-profiles.eps & .png: geographically-arranged subplots of ocean temp & iceberg depth timeseries
%4) Antarctic-iceberg-oceandata-map.eps & .png: map of study site mean iceberg melt rates & ocean temps
%5) *-iceberg-oceandata-scatterplot.eps & .png: site-specific & Antarctic-wide plots of ocean temp above freezing & iceberg melt rate

%%dependencies:
%1) wgs2ps.m
%2) cmocean.m

%THIS README NEEDS TO BE UPDATED

%%%

%% initialize
clearvars;
addpath('/users/ellynenderlin/miscellaneous/general-code','/users/ellynenderlin/miscellaneous/general-code/cmocean');

%specify paths & file names for data
CTD_path = '/Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/CTD_Antarctica/';
CTD_data = [CTD_path,'Antarctic-ocean-data.mat'];
RACMO_path = '/Users/ellynenderlin/Research/miscellaneous/RACMO2.3_Antarctica/';
iceberg_path = '/Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt/';
figure_path = [iceberg_path,'figures/'];

%specify study site names
region = [{'Edgeworth-LarsenA'},{'Crane-LarsenB'},{'Ronne'},{'Filchner'},{'Amery'},{'Totten'},{'Mertz'},...
    {'Thwaites'},{'Ferrigno-Eltanin'},{'Seller-Bugge'},{'Heim-Marguerite'},{'Widdowson-Biscoe'},{'Cadman-Palmer'},{'Blanchard-Danco'},{'Leonardo-Danco'}];
leg_names = [{'Edgeworth'},{'Crane'},{'Ronne'},{'Filchner'},{'Polar Times'},{'Totten'},{'Mertz'},{'Thwaites'},{'Ferrigno'},{'Seller'},{'Heim'},{'Widdowson'},{'Cadman'},{'Blanchard'},{'Leonardo'}];

%specify plot params
marker = ['s','s','s','s','s','s','s','s','s','s','s','s','s','s','s']; %modify if you want to change symbols to indicate something about data (E vs W for example)
plot_letters = [{'i)'},{'j)'},{'k)'},{'l)'},{'m)'},{'n)'},{'o)'},{'h)'},{'g)'},{'f)'},{'e)'},{'d)'},{'c)'},{'b)'},{'a)'}]; %plot letters for sites to be used in geographically-arranged subplots
plot_loc = [2,4,6,8,10,12,14,15,13,11,9,7,5,3,1];
region_colors = [77,172,38; 77,172,38; 184,225,134; 184,225,134; 184,225,134; 184,225,134; 184,225,134;...
    241,182,218; 241,182,218; 208,28,139; 208,28,139; 208,28,139; 208,28,139; 208,28,139; 208,28,139]./255; 
Temp_cmap = cmocean('thermal',600); cmap_add = 3; cmap_mult = 100;
Tforcing_cmap = cmocean('amp',600);

%specify generic variables
rho_sw = 1026; %sea water density in kg m^-3
depth_cutoff = 800; %maximum depth of icebergs & therefore ocean observation data of interest in meters
years = [2011.75 2022.25]; year_ticks = [2013:2:2022]; %approximate date range for plots



%load the RAMP image to plot as background for a map
cd /Users/ellynenderlin/Research/miscellaneous/RAMP
[A,S] = readgeoraster('Antarctic_RAMP_image_v2_1km.tif');
IM.x = S.XWorldLimits(1)+0.5*S.CellExtentInWorldX:S.CellExtentInWorldX:S.XWorldLimits(2)-0.5*S.CellExtentInWorldX;
IM.y = S.YWorldLimits(2)-0.5*S.CellExtentInWorldY:-S.CellExtentInWorldY:S.YWorldLimits(1)+0.5*S.CellExtentInWorldY;
IM.z=single(A);
clear A S;

%navigate to the iceberg directory as the default workspace
cd(iceberg_path);
close all;

%% organize ocean data & pair with remotely-sensed iceberg melt data (run once to create Antarctic-icebergmelt-comparison.mat)
% If you rerun this section, you will overwrite
% Antarctic-icebergmelt-comparison.mat... rerun if you want to change the
% size of the "buffer" search window for ocean data

%specify the buffer region to search for ocean data around the icebergs
buffer = 100000; %m

%load the compiled ocean data
load(CTD_data);

%load RACMO air temps
cd(RACMO_path);
AP_lat = ncread('RACMO2.3p2_XPEN055_T2m_daily_2011_2016.nc','lat');
AP_lon = ncread('RACMO2.3p2_XPEN055_T2m_daily_2011_2016.nc','lon');
AP_airtemp = ncread('RACMO2.3p2_XPEN055_T2m_daily_2011_2016.nc','t2m'); AP_airtemp(AP_airtemp==-9999)=NaN;
AIS_lat = ncread('RACMO2.3p2_ANT27_T2m_daily_2011_2016.nc','lat');
AIS_lon = ncread('RACMO2.3p2_ANT27_T2m_daily_2011_2016.nc','lon');
AIS_airtemp = ncread('RACMO2.3p2_ANT27_T2m_daily_2011_2016.nc','t2m'); AIS_airtemp(AIS_airtemp==-9999)=NaN;
RACMO_time = ncread('RACMO2.3p2_XPEN055_T2m_daily_2011_2016.nc','time');
RACMO_days = RACMO_time-RACMO_time(1); %days since 20110101
for i = 1:size(AP_lat,1); for j = 1:size(AP_lat,2); [AP_x(i,j),AP_y(i,j)] = wgs2ps(AP_lon(i,j),AP_lat(i,j),'StandardParallel',-71,'StandardMeridian',0); end; end
for i = 1:size(AIS_lat,1); for j = 1:size(AIS_lat,2); [AIS_x(i,j),AIS_y(i,j)] = wgs2ps(AIS_lon(i,j),AIS_lat(i,j),'StandardParallel',-71,'StandardMeridian',0); end; end


%loop through the remotely-sensed iceberg melt data & extract data
cd(iceberg_path);
files = dir;
for i = 1:length(region) %loop through the regions
    %find the folder containing the data for the specified region
    for j = 1:length(files)
        if ~isempty(strfind(files(j).name,string(region(i)))) && isdir(files(j).name)
            regionref = j;
        end
    end
    
    % extract the data for the region
    cd_to_region = ['cd ',files(regionref).name]; eval(cd_to_region);
    melt(i).name = region(i); melt(i).dispname = leg_names(i);
    meltinfo = dir('*iceberg_meltinfo.csv');
    date_o = []; xcoord_o = []; ycoord_o = []; 
    date_f = []; xcoord_f = []; ycoord_f = [];
    rhoi_o = []; rhoi_f = [];
    flux = []; meltrate = []; meltrate_sigma = []; 
    keeld = []; keeld_sigma = []; 
    surf_area = []; surf_area_sigma = []; 
    sub_area = []; sub_area_sigma = []; 
    for k = 1:length(meltinfo)
        M=readtable(meltinfo(k).name); %table exported using plot_export_iceberg_melt_data.m
        M = table2array(M); %convert to array to enable easier indexing
        
        %pull variables
        xo = M(:,2); yo = M(:,3); zo = M(:,4); rhoo = M(:,5); %initial locations, median elev, density
        xf = M(:,7); yf = M(:,8); zf = M(:,9); rhof = M(:,10); %same as above but final
        dVdt = M(:,16); dVdt_uncert = M(:,17); %m^3/d

        %recalculate draft and Asub
        draft = (nanmean([rhoo rhof],2)./(repmat(rho_sw,length(nanmean([rhoo rhof],2)),1)-nanmean([rhoo rhof],2))).*nanmean([zo zf],2); %draft = M(:,18); 
        draft_uncert = M(:,19);
        Asurf = M(:,20); Asurf_uncert = M(:,21);
        lat_area = M(:,22) - Asurf; perim = lat_area./draft; clear lat_area; lat_area = perim.*draft; Asub = lat_area + Asurf; clear lat_area perim; %Asub = M(:,22); 
        Asub_uncert = M(:,23);
        rhoi_o = [rhoi_o; rhoo]; rhoi_f = [rhoi_f; rhof];
        surf_area = [surf_area; Asurf]; surf_area_sigma = [surf_area_sigma; Asurf_uncert]; 
        sub_area = [sub_area; Asub]; sub_area_sigma = [sub_area_sigma; Asub_uncert]; 
        meltrate = [meltrate; (dVdt./Asub)]; meltrate_sigma = [meltrate_sigma; abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2)]; 
        keeld = [keeld; draft]; keeld_sigma = [keeld_sigma; draft_uncert];
        
        %converting to decimal dates
        decidate = convert_to_decimaldate(meltinfo(k).name(end-49:end-42));
        date_o = [date_o; repmat(decidate,size(xo))]; xcoord_o = [xcoord_o; xo]; ycoord_o = [ycoord_o; yo]; clear decidate;
        decidate = convert_to_decimaldate(meltinfo(k).name(end-34:end-27));
        date_f = [date_f; repmat(decidate,size(xf))]; xcoord_f = [xcoord_f; xf]; ycoord_f = [ycoord_f; yf]; clear decidate;
        
        clear xo yo zo rhoo xf yf zf rhof dVdt* draft* Asurf* Asub* M;
    end
    melt(i).to = date_o; melt(i).tf = date_f;
    melt(i).x = nanmean([xcoord_o xcoord_f],2); melt(i).y = nanmean([ycoord_o ycoord_f],2); 
    melt(i).rho = nanmean([rhoi_o rhoi_f],2); 
    melt(i).d = keeld; melt(i).d_uncert = keeld_sigma; 
    melt(i).m = meltrate; melt(i).m_uncert = meltrate_sigma; 
    melt(i).Asurf = surf_area; melt(i).Asurf_uncert = surf_area_sigma; 
    melt(i).Asub = sub_area; melt(i).Asub_uncert = sub_area_sigma; 
    avgx(i) = nanmean(melt(i).x); avgy(i) = nanmean(melt(i).y); 
    clear date_* *coord_* keeld meltrate *_uncert;
    cd .. %back out of the region directory
    
    %extract RACMO-estimated air temperatures from the nearest grid cell &
    %average over the full record to estimate iceberg temperature
    if i <= 2 || i >= 10
        lat_diff = abs(repmat(nanmean(melt(i).y),size(AP_y)) - AP_y);
        lon_diff = abs(repmat(nanmean(melt(i).x),size(AP_x)) - AP_x);
        diff_map = sqrt(lat_diff.^2+lon_diff.^2); diff_map(isnan(squeeze(nanmean(AP_airtemp(:,:,1,270:450),4)))) = NaN; %solve for the distance vector using the x&y distances
        RACMO_ref = find(diff_map==min(min(diff_map))); %find the minimum distance (reference for your grid cell is output)
        [RACMOy RACMOx] = ind2sub(size(squeeze(nanmean(AP_airtemp(:,:,1,270:450),4))),RACMO_ref); %convert cell reference to an x- and y-cell index
        melt(i).bergT = squeeze(nanmean(AP_airtemp(RACMOx,RACMOy,1,:),4))-273.15;
    else
        lat_diff = abs(repmat(nanmean(melt(i).y),size(AIS_y)) - AIS_y);
        lon_diff = abs(repmat(nanmean(melt(i).x),size(AIS_x)) - AIS_x);
        diff_map = sqrt(lat_diff.^2+lon_diff.^2); diff_map(isnan(squeeze(nanmean(AIS_airtemp(:,:,1,270:450),4)))) = NaN; %solve for the distance vector using the x&y distances
        RACMO_ref = find(diff_map==min(min(diff_map))); %find the minimum distance (reference for your grid cell is output)
        [RACMOy RACMOx] = ind2sub(size(squeeze(nanmean(AIS_airtemp(:,:,1,270:450),4))),RACMO_ref); %convert cell reference to an x- and y-cell index
        melt(i).bergT = squeeze(nanmean(AIS_airtemp(RACMOx,RACMOy,1,:),4))-273.15;
    end
    clear *diff* RACMOy RACMOx;
    
    %crop the CTD and APB data to within the specified buffer distance of the icebergs
    xlimits = [min(melt(i).x)-buffer max(melt(i).x)+buffer]; ylimits = [min(melt(i).y)-buffer max(melt(i).y)+buffer];
    %CTD
    for j = 1:length(CTD)
        %find data
        oceanrefs = find(CTD(j).x>=min(xlimits) & CTD(j).x<=max(xlimits) & CTD(j).y>=min(ylimits) & CTD(j).y<=max(ylimits));
        if ~isempty(oceanrefs)
            melt(i).oceant = CTD(j).date(oceanrefs);
            melt(i).oceanx = CTD(j).x(oceanrefs); melt(i).oceany = CTD(j).y(oceanrefs);
            melt(i).oceand = CTD(j).depth(:,oceanrefs); melt(i).oceanT = CTD(j).T(:,oceanrefs); melt(i).oceanS = CTD(j).sal(:,oceanrefs);
            
            %if along peninsula, filter across the spine
            if avgx(i) < -1.8e6
                if avgx(i)>-2.45e6
                    melt(i).oceant(melt(i).oceanx<-2.45e6) = [];
                    melt(i).oceand(:,melt(i).oceanx<-2.45e6) = [];
                    melt(i).oceanT(:,melt(i).oceanx<-2.45e6) = [];
                    melt(i).oceanS(:,melt(i).oceanx<-2.45e6) = [];
                    melt(i).oceany(melt(i).oceanx<-2.45e6) = [];
                    melt(i).oceanx(melt(i).oceanx<-2.45e6) = [];
                else
                    melt(i).oceant(melt(i).oceanx>-2.45e6) = [];
                    melt(i).oceand(:,melt(i).oceanx>-2.45e6) = [];
                    melt(i).oceanT(:,melt(i).oceanx>-2.45e6) = [];
                    melt(i).oceanS(:,melt(i).oceanx>-2.45e6) = [];
                    melt(i).oceany(melt(i).oceanx>-2.45e6) = [];
                    melt(i).oceanx(melt(i).oceanx>-2.45e6) = [];
                end
            end
            
        end
        if ~isempty(oceanrefs)
            disp(['Number of CTD observations for ',char(region(i)),' = ',num2str(length(oceanrefs))]);
        end
        clear oceanrefs;
    end
    %APB
    for j = 1:length(APB)
        %find data
        oceanrefs = find(APB(j).x>=min(xlimits) & APB(j).x<=max(xlimits) & APB(j).y>=min(ylimits) & APB(j).y<=max(ylimits));
        if ~isempty(oceanrefs)
            melt(i).oceant = APB(j).date(oceanrefs);
            melt(i).oceanx = APB(j).x(oceanrefs); melt(i).oceany = APB(j).y(oceanrefs);
            melt(i).oceand = APB(j).depth(:,oceanrefs); melt(i).oceanT = APB(j).T(:,oceanrefs); melt(i).oceanS = APB(j).sal(:,oceanrefs);
            
            %if along peninsula, filter across the spine
            if avgx(i) < -1.8e6
                if avgx(i)>-2.45e6
                    melt(i).oceant(melt(i).oceanx<-2.45e6) = [];
                    melt(i).oceand(:,melt(i).oceanx<-2.45e6) = [];
                    melt(i).oceanT(:,melt(i).oceanx<-2.45e6) = [];
                    melt(i).oceanS(:,melt(i).oceanx<-2.45e6) = [];
                    melt(i).oceany(melt(i).oceanx<-2.45e6) = [];
                    melt(i).oceanx(melt(i).oceanx<-2.45e6) = [];
                else
                    melt(i).oceant(melt(i).oceanx>-2.45e6) = [];
                    melt(i).oceand(:,melt(i).oceanx>-2.45e6) = [];
                    melt(i).oceanT(:,melt(i).oceanx>-2.45e6) = [];
                    melt(i).oceanS(:,melt(i).oceanx>-2.45e6) = [];
                    melt(i).oceany(melt(i).oceanx>-2.45e6) = [];
                    melt(i).oceanx(melt(i).oceanx>-2.45e6) = [];
                end
            end
            
        end
        if ~isempty(oceanrefs)
            disp(['Number of APB observations for ',char(region(i)),' = ',num2str(length(oceanrefs))]);
        end
        clear oceanrefs;
    end
    clear *limits;
    
end
save([iceberg_path,'Antarctic-icebergmelt-comparison.mat'],'melt','-v7.3');
clear CTD APB *airtemp;
disp('Done compiling iceberg melt estimates, RACMO air temps (for iceberg temps), & ocean data');


% determine the appropriate depth resolution of the ocean observations
close all; warning off;
figure; set(gcf,'position',[50 50 800 400]);
sub1 = subplot(1,2,1); sub2 = subplot(1,2,2);
for i = 1:length(melt)
    disp(char(melt(i).name));
    
    %calculate the median profile depth increment for each profile
    if size(melt(i).oceand,2) > 0
        for j = 1:size(melt(i).oceand,2)
            dz0_200m(j) = nanmedian(diff(melt(i).oceand(find(melt(i).oceand(:,j)<=200),j)));
            dz201_800m(j) = nanmedian(diff(melt(i).oceand(find(melt(i).oceand(:,j)>200 & melt(i).oceand(:,j)<=depth_cutoff),j)));
            dz0_800m(j) = nanmedian(diff(melt(i).oceand(find(melt(i).oceand(:,j)<=depth_cutoff),j)));
        end
        
        %plot histograms of median profile increments
        subplot(sub1); hps(i) = histogram(dz0_200m); hold on;
        subplot(sub2); hpd(i) = histogram(dz201_800m); hold on;
        
        %print the median profile increment for the site
        fprintf('median depth increment for shallow water: %i \n',round(nanmedian(dz0_200m)));
        fprintf('median depth increment for deep water: %i \n',round(nanmedian(dz201_800m)));
        
        %create a standard depth profile for the site
        melt(i).oceandavg_prof = [0:round(nanmedian(dz0_800m)):depth_cutoff]';
        
        clear dz0_200m dz201_800m dz0_800m;
    else
        melt(i).oceandavg_prof = [0:10:depth_cutoff]';
        fprintf('no ocean observations for site \n');
    end
    subplot(sub1); ylabel('Count'); xlabel('Binning increment for depths <200 m (m)');
    subplot(sub2); xlabel('Binning increment for depths 201-800 m (m)');
end
drawnow;
save([iceberg_path,'Antarctic-icebergmelt-comparison.mat'],'melt','-v7.3');

%% interpolate profiles to a standard, uniform, site-specific profile
warning off;
if ~exist('melt')
    load([iceberg_path,'Antarctic-icebergmelt-comparison.mat']);
end
% melt = rmfield(melt,{'oceand_prof','oceanTavg_prof','oceanTstd_prof','oceanSavg_prof','oceanSstd_prof'}); melt = rmfield(melt,'bergtemp');

%interpolate data to a standard, regular depth profile and then temporally
%average if multiple profiles were collected for the same date
disp('Interpolating to a standard uniform depth scale based on median depth increment of data for site');
for i = 1:length(melt)
    if ~isempty(melt(i).oceant)
        disp(char(melt(i).name));
        
        %create average temperature & salinity profiles for each date
        [unique_dates,unique_refs] = unique(melt(i).oceant);
        for j = 1:length(unique_dates)
            year_refs = find(melt(i).oceant == unique_dates(j));
            melt(i).oceantavg(j) = unique_dates(j);
            if length(year_refs) > 1
                melt(i).oceanxavg(:,j) = nanmean(melt(i).oceanx(year_refs)); melt(i).oceanyavg(:,j) = nanmean(melt(i).oceany(year_refs));
                for k = 1:length(year_refs)
                    if sum(~isnan(melt(i).oceand(:,year_refs(k))))>1 && sum(~isnan(melt(i).oceanT(:,year_refs(k))))>1
                        Tprof(:,k) = interp1(melt(i).oceand(~isnan(melt(i).oceand(:,year_refs(k))) & ~isnan(melt(i).oceanT(:,year_refs(k))),year_refs(k)),melt(i).oceanT(~isnan(melt(i).oceand(:,year_refs(k))) & ~isnan(melt(i).oceanT(:,year_refs(k))),year_refs(k)),melt(i).oceandavg_prof,'linear');
                        Sprof(:,k) = interp1(melt(i).oceand(~isnan(melt(i).oceand(:,year_refs(k))) & ~isnan(melt(i).oceanT(:,year_refs(k))),year_refs(k)),melt(i).oceanS(~isnan(melt(i).oceand(:,year_refs(k))) & ~isnan(melt(i).oceanT(:,year_refs(k))),year_refs(k)),melt(i).oceandavg_prof,'linear');
                    else
                        Tprof(:,k) = NaN(size(melt(i).oceandavg_prof)); Sprof(:,k) = NaN(size(melt(i).oceandavg_prof));
                    end
                end
                melt(i).oceanTavg_prof(:,j) = nanmean(Tprof,2); melt(i).oceanTstd_prof(:,j) = nanstd(Tprof,0,2);
                melt(i).oceanSavg_prof(:,j) = nanmean(Sprof,2); melt(i).oceanSstd_prof(:,j) = nanstd(Sprof,0,2);
                melt(i).oceanTfreeze_prof(:,j) = (((-5.73*10^-2).*melt(i).oceanSavg_prof(:,j)) + (8.32*10^-2) - ((7.61*10^-4).*melt(i).oceandavg_prof));
                clear Tprof Sprof;
            else
                melt(i).oceanxavg(:,j) = melt(i).oceanx(year_refs); melt(i).oceanyavg(:,j) = melt(i).oceany(year_refs);
                if sum(~isnan(melt(i).oceand(:,year_refs)))>1 && sum(~isnan(melt(i).oceanT(:,year_refs)))>1
                    melt(i).oceanTavg_prof(:,j) = interp1(melt(i).oceand(~isnan(melt(i).oceand(:,year_refs)) & ~isnan(melt(i).oceanT(:,year_refs)),year_refs),melt(i).oceanT(~isnan(melt(i).oceand(:,year_refs)) & ~isnan(melt(i).oceanT(:,year_refs)),year_refs),melt(i).oceandavg_prof,'linear'); 
                    melt(i).oceanSavg_prof(:,j) = interp1(melt(i).oceand(~isnan(melt(i).oceand(:,year_refs)) & ~isnan(melt(i).oceanT(:,year_refs)),year_refs),melt(i).oceanS(~isnan(melt(i).oceand(:,year_refs)) & ~isnan(melt(i).oceanT(:,year_refs)),year_refs),melt(i).oceandavg_prof,'linear'); 
                else
                    melt(i).oceanTavg_prof(:,j) = NaN(size(melt(i).oceandavg_prof)); melt(i).oceanSavg_prof(:,j) = NaN(size(melt(i).oceandavg_prof));
                end
                melt(i).oceanTstd_prof(:,j) = NaN(size(melt(i).oceandavg_prof)); melt(i).oceanSstd_prof(:,j) = NaN(size(melt(i).oceandavg_prof));
                melt(i).oceanTfreeze_prof(:,j) = (((-5.73*10^-2).*melt(i).oceanSavg_prof(:,j)) + (8.32*10^-2) - ((7.61*10^-4).*melt(i).oceandavg_prof));
            end

            clear year_refs;
        end
        clear unique_*;
        
        %add date-averaged profiles to structure
        save([iceberg_path,'Antarctic-icebergmelt-comparison.mat'],'melt','-v7.3');
    end
    
end


%% create profiles of ocean temperature & iceberg depth timeseries (Antarctic-iceberg-oceandata-profiles.eps)

%reload compiled data as needed
if ~exist('melt')
    load([iceberg_path,'Antarctic-icebergmelt-comparison.mat']);
end

%set up plot
close all; 
figure; set(gcf,'position',[50 50 800 1000]);
imagesc(IM.x,IM.y,IM.z); colormap gray; hold on; axis xy equal; colormap(gray(10001));

%loop through data
disp('Plotting temperature profiles & iceberg depth timeseries');
for i = 1:length(melt)
    if ~isempty(melt(i).oceant)
        subpl = subplot(8,2,plot_loc(i));
        
        %plot temperature profiles with colors to distinguish temps
        for j = 1:length(melt(i).oceantavg)
            prof_cmap = []; prof_d = [];
            
            if nanmedian(diff(melt(i).oceandavg_prof)) == 1
                %average over 20 rows to create smoothed profiles at 20 m increments
                for k = 11:20:size(melt(i).oceanTavg_prof(:,j),1)-10
                    prof_d = [prof_d; nanmean(melt(i).oceandavg_prof(k-10:k+10))];
                    if ~isnan(nanmean(melt(i).oceanTavg_prof(k-10:k+10,j))) %&& ~isnan(melt(i).oceanTfreeze_prof(k,j))
                        prof_cmap = [prof_cmap; Temp_cmap(round((nanmean(melt(i).oceanTavg_prof(k-10:k+10,j))+cmap_add)*cmap_mult),:)];
                        %                         plot(melt(i).oceantavg(j),nanmean(melt(i).oceandavg_prof(k-10:k+10)),'.','color',Temp_cmap(round((nanmean(melt(i).oceanTavg_prof(k-10:k+10,j))+3)*100),:),'markersize',10); hold on;
                    else
                        prof_cmap = [prof_cmap; 1,1,1];
                    end
                end
            else
                %plot at native standardized resolution because it should be ~20 m
                for k = 1:size(melt(i).oceanTavg_prof(:,j),1)
                    prof_d = [prof_d; melt(i).oceandavg_prof(k)];
                    if ~isnan(melt(i).oceanTavg_prof(k,j))
                        prof_cmap = [prof_cmap; Temp_cmap(round((melt(i).oceanTavg_prof(k,j)+cmap_add)*cmap_mult),:)];
                        %                         plot(melt(i).oceantavg(j),melt(i).oceandavg_prof(k),'.','color',Temp_cmap(round((melt(i).oceanTavg_prof(k,j)+3)*100),:),'markersize',10); hold on;
                    else
                        prof_cmap = [prof_cmap; 1,1,1];
                    end
                end
            end
            top_ref = find(sum(prof_cmap,2) < 3,1,'first'); %identify the deepest part of the profile with data as the addition of the plot colors < 3 (sum = 3 is white for NaNs)
            bottom_ref = find(sum(prof_cmap,2) < 3,1,'last'); %identify the deepest part of the profile with data as the addition of the plot colors < 3 (sum = 3 is white for NaNs)
            scatter(repmat(melt(i).oceantavg(j),size(prof_d(top_ref:bottom_ref))),prof_d(top_ref:bottom_ref),14,prof_cmap(top_ref:bottom_ref,:),'filled','s'); hold on;
            clear prof_cmap prof_d bottom_ref;
        end
        
        %plot the iceberg depths for each date
        plot(nanmean([melt(i).to melt(i).tf],2),melt(i).d,'.k'); hold on;
        errorbar(nanmean([melt(i).to melt(i).tf],2),melt(i).d,[],[],abs(nanmean([melt(i).to melt(i).tf],2)-melt(i).to),abs(nanmean([melt(i).to melt(i).tf],2)-melt(i).tf),'.k');
        
        %format plot
        set(gca,'ydir','reverse','xlim',years,'ylim',[0 depth_cutoff],'ytick',[0:250:750],'fontsize',16);
        if plot_loc(i) == 15
            xlabel('Year','fontsize',16); ylabel('Depth (m b.s.l.)','fontsize',16);
        end
        text(min(get(gca,'xlim'))+0.025*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),max(get(gca,'ylim'))-0.2*abs(max(get(gca,'ylim'))-min(get(gca,'ylim'))),[char(plot_letters(i)),' ',char(melt(i).dispname)],'fontsize',16);
        drawnow;

    else
        subpl = subplot(8,2,plot_loc(i));
        %plot the iceberg depths for each date
        plot(nanmean([melt(i).to melt(i).tf],2),melt(i).d,'.k'); hold on;
        errorbar(nanmean([melt(i).to melt(i).tf],2),melt(i).d,[],[],abs(nanmean([melt(i).to melt(i).tf],2)-melt(i).to),abs(nanmean([melt(i).to melt(i).tf],2)-melt(i).tf),'.k');
        
        %format plot
        set(gca,'ydir','reverse','xlim',years,'ylim',[0 depth_cutoff],'ytick',[0:250:750],'fontsize',16);
        text(min(get(gca,'xlim'))+0.025*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),max(get(gca,'ylim'))-0.2*abs(max(get(gca,'ylim'))-min(get(gca,'ylim'))),[char(plot_letters(i)),' ',char(melt(i).dispname)],'fontsize',16);
        drawnow;
        
    end
    
    %format axes
    pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.05*pos(3) 1.15*pos(4)]);
    if plot_loc(i) == 14
        set(gca,'xlim',years,'xtick',year_ticks,'ylim',[0 depth_cutoff],'ytick',[0:250:750],'fontsize',16);
        xlabel('Year','fontsize',16); 
    elseif plot_loc(i) == 15
        set(gca,'xlim',years,'xtick',year_ticks,'ylim',[0 depth_cutoff],'ytick',[0:250:750],'fontsize',16);
        xlabel('Year','fontsize',16); ylbl = ylabel('Depth (m b.s.l.)','fontsize',16);
        set(ylbl,'position',[min(years)-1.5 -3000 -1]);
    else
        set(gca,'xlim',years,'xtick',year_ticks,'xticklabel',[],'ylim',[0 depth_cutoff],'ytick',[0:250:750],'fontsize',16);
    end
    box on;
    
end
%add a colorbar
annotation('rectangle',[0.57 0.11 0.35 0.065],'facecolor','w','edgecolor','k');
for j = 1:length(Temp_cmap)
    annotation('line',[0.595+j/2000 0.595+j/2000],[0.1525 0.170],'color',Temp_cmap(j,:),'linewidth',1.5);
end
annotation('textbox',[0.57 0.135 0.05 0.02],'string',['-',num2str(cmap_add),char(176),'C'],'fontsize',16,'edgecolor','none');
annotation('textbox',[0.725 0.135 0.05 0.02],'string',['0',char(176),'C'],'fontsize',16,'edgecolor','none');
annotation('textbox',[0.875 0.135 0.05 0.02],'string',[num2str(cmap_add),char(176),'C'],'fontsize',16,'edgecolor','none');
annotation('textbox',[0.65 0.115 0.25 0.02],'string','ocean temperature','fontsize',16,'edgecolor','none');

%save
saveas(gcf,[figure_path,'Antarctic-iceberg-oceandata-profiles.eps'],'epsc'); saveas(gcf,[figure_path,'Antarctic-iceberg-oceandata-profiles.png'],'png');
disp('iceberg and ocean temp depth profiles saved');

% %create a map that shows the average temperature for each profile down to
% %~100 m depth (Xs) and the median iceberg depth for all sites
% figure(Tm_mapplot);
% for i = 1:length(melt)
%     if ~isempty(melt(i).oceant)
%         for j = 1:length(melt(i).oceantavg)
%             hundred_ref = find(melt(i).oceandavg_prof<=100,1,'last');
%             median_ref = find(melt(i).oceandavg_prof<=nanmedian(melt(i).d),1,'last');
%             if ~isnan(nanmean(melt(i).oceanTavg_prof(1:hundred_ref,j)-melt(i).oceanTfreeze_prof(1:hundred_ref,j)))
%             plot(melt(i).oceanxavg(j),melt(i).oceanyavg(j),'x','color',Temp_cmap(round((nanmean(melt(i).oceanTavg_prof(1:hundred_ref,j)-melt(i).oceanTfreeze_prof(1:hundred_ref,j))+1)*100),:)); hold on;
%             end
%         end
%     end
%     plot(nanmean(melt(i).x),nanmean(melt(i).y),[marker(i),'k'],'markerfacecolor',depth_cmap(round(nanmean(melt(i).d)),:),'markersize',12); hold on;
% end
% %add labels to the location plot
% for i = 1:length(melt)
%     figure(Tm_mapplot); 
%     if strcmp(marker(i),'d')
%         text(nanmean(melt(i).x)+100000,nanmean(melt(i).y),char(plot_letters(i)),'fontsize',16);
%     else
%         if strcmp(char(plot_letters(i)),'f)')
%             text(nanmean(melt(i).x)-200000,nanmean(melt(i).y)-100000,char(plot_letters(i)),'fontsize',16); 
%         elseif strcmp(char(plot_letters(i)),'a)')
%             text(nanmean(melt(i).x)-200000,nanmean(melt(i).y)+100000,char(plot_letters(i)),'fontsize',16); 
%         else
%             text(nanmean(melt(i).x)-200000,nanmean(melt(i).y),char(plot_letters(i)),'fontsize',16); 
%         end
%     end
% end
% %label
% set(gca,'xlim',[-28e5 28e5],'xtick',[-24e5:8e5:24e5],'xticklabel',[-2400:800:2400],...
%     'ylim',[-24e5 24e5],'ytick',[-24e5:8e5:24e5],'yticklabel',[-2400:800:2400],'fontsize',24); grid off;
% xlabel('Easting (km)','fontsize',24); ylabel('Northing (km)','fontsize',24);
% graticuleps(-50:-5:-90,-180:30:180);
% text(0,6.5e5,'85^oS','fontsize',16); text(0,12.0e5,'80^oS','fontsize',16); text(0,17.5e5,'75^oS','fontsize',16); text(0,23.0e5,'70^oS','fontsize',16);
% text(-16.5e5,25.25e5,'-30^oE','fontsize',16); text(12.5e5,25.25e5,'30^oE','fontsize',16); 
% colormap(gca,gray(100001));
% saveas(gcf,'Antarctic-iceberg-oceandata-map.eps','epsc'); saveas(gcf,'Antarctic-iceberg-oceandata-map.png','png');
% %now zoom in on the peninsula and save again
% set(gca,'xlim',[-28e5 -20e5],'xtick',[-28e5:2e5:-20e5],'xticklabel',[-2800:200:-2000],...
%     'ylim',[7.5e5 17.5e5],'ytick',[8e5:2e5:16e5],'yticklabel',[800:200:1600],'fontsize',24);
% graticuleps(-50:-2:-90,-180:10:180);
% saveas(gcf,'AntarcticPeninsula-iceberg-oceandata-map.eps','epsc'); saveas(gcf,'AntarcticPeninsula-iceberg-oceandata-map.png','png');

%% compute thermal forcing & generate plots to show melt rate and thermal forcing together
close all;
disp('Solving for thermal forcing & plotting figures to show iceberg melt rates vs thermal forcing');

%reload compiled data as needed
if ~exist('melt')
    load([iceberg_path,'Antarctic-icebergmelt-comparison.mat']);
end

%set-up a standard depth profile
Tm_scatterplot = figure; set(gcf,'position',[850 50 800 400]); 
depth_cmap = cmocean('deep',depth_cutoff); im_cmap = colormap(gray(10001)); im_cmap(1,:) = [1 1 1];
Tm_mapplot = figure; set(gcf,'position',[50 50 800 800]);
imagesc(IM.x,IM.y,IM.z); hold on; axis xy equal; colormap(gca,im_cmap);
tempref = [];
for i = 1:size(region,2)
    disp_names(i) = {strjoin([cellstr(plot_letters(i)),' ',cellstr(leg_names(i))])};
    landsats = dir([iceberg_path,char(region(i)),'/LC*PS.TIF']);
    sitex = []; sitey = []; %set up empty cells to insert ocean data coordinates (if data exist) & iceberg coordinates for adjusting site maps
    
    %map for study site
    Tmsitemap = figure; set(gcf,'position',[50 850 800 450]);
    [A,S] = readgeoraster([landsats(1).folder,'/',landsats(1).name]);
    im.x = S.XWorldLimits(1):S.SampleSpacingInWorldX:S.XWorldLimits(2);
    im.y = S.YWorldLimits(2):-S.SampleSpacingInWorldY:S.YWorldLimits(1);
    im.z = double(A); clear A S landsats;
    imagesc(im.x,im.y,im.z); axis xy equal; colormap(gray(10001)); hold on;
    
    %find appropriate ocean data & extract thermal forcing estimates
    if ~isempty(melt(i).oceant)
        TFcoords = []; TF = [];
        
        %loop through remotely-sensed data and extract ocean forcing information for each iceberg
        for j = 1:length(melt(i).to) %length of melt(i).to corresponds to the number of icebergs
            %identify the time span of remotely-sensed iceberg melt rate estimates
            %(bi-annual=2, annual=1, or seasonal=0)
            if melt(i).tf(j)-melt(i).to(j) >= 2
                timespan = 2;
            elseif melt(i).tf(j)-melt(i).to(j) >= 1
                timespan = 1;
            else
                timespan = 0;
            end
            
            %if seasonal, find ocean data from approximately the same season
            %minrefs specifies the indices for the closest date (there may be multiple profiles) & oceantemps and oceansals are the corresponding profiles
            if size(melt(i).oceanx,2) ~=1; melt(i).oceanx = melt(i).oceanx'; melt(i).oceany = melt(i).oceany'; end %make sure coordinates are always a column vector
            if timespan == 0
                deciseas = nanmean([melt(i).to(j) melt(i).tf(j)]-floor(melt(i).to(j)),2); if deciseas > 1; deciseas = deciseas - floor(deciseas); end
                
                [mindiff,minref] = min(abs((melt(i).oceant-floor(melt(i).oceant))-deciseas));
                if melt(i).oceant(minref)-floor(melt(i).oceant(minref)) >= melt(i).to(j)-floor(melt(i).to(j)) && melt(i).oceant(minref)-floor(melt(i).oceant(minref)) <= melt(i).tf(j)-floor(melt(i).to(j)) %if to and tf are in the same year & minref is in between, find all between
                    minrefs = find(melt(i).oceant-floor(melt(i).oceant(minref)) >= melt(i).to(j)-floor(melt(i).to(j)) & melt(i).oceant-floor(melt(i).oceant(minref)) <= melt(i).tf(j)-floor(melt(i).to(j)));
                elseif melt(i).oceant(minref)-floor(melt(i).oceant(minref)) <= melt(i).to(j)-floor(melt(i).to(j)) && melt(i).oceant(minref)-floor(melt(i).oceant(minref)) <= melt(i).tf(j)-floor(melt(i).tf(j)) %if tf is in a different year than to & minref is in between, find all between
                    minrefs = find(melt(i).oceant-floor(melt(i).oceant(minref)) <= melt(i).to(j)-floor(melt(i).to(j)) & melt(i).oceant-floor(melt(i).oceant(minref)) <= melt(i).tf(j)-floor(melt(i).tf(j)));
                else
                    if mindiff < 0.5 %if there are no data that fall within the seasonal range of to and tf, find data within +/-3 months of the central day of year
                        minrefs = find(abs((melt(i).oceant-floor(melt(i).oceant))-deciseas) <= 4/12);
                    else
                        minrefs = find(abs((melt(i).oceant-floor(melt(i).oceant))-deciseas) <= mindiff + 1/12);
                    end
                end
                melt(i).oceanm_tdiff(j,:) = [nanmean([melt(i).to(j) melt(i).tf(j)])-min(melt(i).oceant(minrefs)) max(melt(i).oceant(minrefs))-nanmean([melt(i).to(j) melt(i).tf(j)])];
                oceanx = melt(i).oceanx(minrefs); oceany = melt(i).oceany(minrefs);
                oceantemps = melt(i).oceanT(:,minrefs); oceansals = melt(i).oceanS(:,minrefs); oceandepths = melt(i).oceand(:,minrefs);
                clear minref deciseas mindiff;
            else
                %if annual or bi-annual, find the closest year of ocean data
                [~,minref] = min(abs(melt(i).oceant-nanmean([melt(i).to(j) melt(i).tf(j)])));
                if melt(i).oceant(minref)-nanmean([melt(i).to(j) melt(i).tf(j)]) > 0
                    minrefs = find(melt(i).oceant>=melt(i).oceant(minref) & melt(i).oceant<=melt(i).oceant(minref)+timespan);
                else
                    minrefs = find(melt(i).oceant<=melt(i).oceant(minref) & melt(i).oceant>=melt(i).oceant(minref)-timespan);
                end
                %calculate the time difference between remotely-sensed and in situ observations
                if min(melt(i).oceant(minrefs)) > melt(i).tf(j)
                    melt(i).oceanm_tdiff(j,:) = [min(melt(i).oceant(minrefs))-nanmean([melt(i).to(j) melt(i).tf(j)]) max(melt(i).oceant(minrefs))-nanmean([melt(i).to(j) melt(i).tf(j)])];
                elseif max(melt(i).oceant(minrefs)) < melt(i).to(j)
                    melt(i).oceanm_tdiff(j,:) = [nanmean([melt(i).to(j) melt(i).tf(j)])-min(melt(i).oceant(minrefs)) nanmean([melt(i).to(j) melt(i).tf(j)])-max(melt(i).oceant(minrefs))];
                else
                    melt(i).oceanm_tdiff(j,:) = [nanmean([melt(i).to(j) melt(i).tf(j)])-min(melt(i).oceant(minrefs)) max(melt(i).oceant(minrefs))-nanmean([melt(i).to(j) melt(i).tf(j)])];
                end
                oceanx = melt(i).oceanx(minrefs); oceany = melt(i).oceany(minrefs);
                oceantemps = melt(i).oceanT(:,minrefs); oceansals = melt(i).oceanS(:,minrefs); oceandepths = melt(i).oceand(:,minrefs);
                clear minref;
                
            end
            
            %extract temperature metrics over the iceberg draft
            for k = 1:length(minrefs)
                %identify the bottom of each profile
                if ~isempty(find(oceandepths(:,k)<=melt(i).d(j),1,'last'))
                    bottomrefs(k) = find(oceandepths(:,k)<=melt(i).d(j),1,'last'); %index for the deepest observation for each profile
                    bottomT(k) = melt(i).oceanT(bottomrefs(k),minrefs(k)); bottomS(k) = melt(i).oceanS(bottomrefs(k),minrefs(k));
                else
                    bottomrefs(k) = NaN; bottomT(k) = NaN; bottomS(k) = NaN;
                end
                
                %use the trapz function to calculate the mean for each profile
                %if its maximum observation depth is >90% of the iceberg draft
                if 0.9*max(oceandepths(~isnan(oceantemps(:,k)),k)) > melt(i).d(j) & min(oceandepths(~isnan(oceantemps(:,k)),k)) < 50
                    Tavg(k) = vertmean2(-oceandepths(:,k),oceantemps(:,k),-melt(i).d(j));
                else
                    Tavg(k) = NaN;
                end
                %repeat averaging but for salinity (may not have salinity corresponding to all temp observations)
                if 0.9*max(oceandepths(~isnan(oceansals(:,k)),k)) > melt(i).d(j) & min(oceandepths(~isnan(oceansals(:,k)),k)) < 50
                    Savg(k) = vertmean2(-oceandepths(:,k),oceansals(:,k),-melt(i).d(j));
                else
                    Savg(k) = NaN;
                end
                %repeat averaging, but to calculate the average freezing temperature of sea water (Tfp = -5.73*10^-2 (C/psu)*salinity + 8.32*10^-2 (C) - 7.61*10^-4 (C/dbar)*pressure)
                %pressure is approximately equivalent to depth
                if 0.9*max(oceandepths(~isnan(oceansals(:,k)),k)) > melt(i).d(j) & min(oceandepths(~isnan(oceansals(:,k)),k)) < 50
                    Tfpavg(k) = vertmean2(-oceandepths(:,k),(((-5.73*10^-2).*oceansals(:,k)) + (8.32*10^-2) - ((7.61*10^-4).*oceandepths(:,k))),-melt(i).d(j));
                else
                    Tfpavg(k) = NaN;
                end
            end
            
            %extract the average, standard deviation, max, & bottom temp and salinity from the profiles for the closest data
            %             melt(i).oceanTavg(j,1) = nanmean(oceantemps(oceandepths<=melt(i).d(j)));
            melt(i).oceanTavg(j,1) = nanmean(Tavg);
            melt(i).oceanTstd(j,1) = nanstd(oceantemps(oceandepths<=melt(i).d(j)));
            melt(i).oceanTmax(j,1) = max(oceantemps(oceandepths<=melt(i).d(j))); melt(i).oceanTbottom(j,1) = nanmean(bottomT);
            %             melt(i).oceanSavg(j,1) = nanmean(oceansals(oceandepths<=melt(i).d(j)));
            melt(i).oceanSavg(j,1) = nanmean(Savg);
            melt(i).oceanSstd(j,1) = nanstd(oceansals(oceandepths<=melt(i).d(j)));
            melt(i).oceanSmax(j,1) = max(oceansals(oceandepths<=melt(i).d(j))); melt(i).oceanSbottom(j,1) = nanmean(bottomS);
            melt(i).oceanTfp(j,1) = nanmean(Tfpavg);
            clear bottom*;
            
            %set the size of the melt rate vs thermal forcing scatterplot
            %symbols so that they vary with draft
            draft_size(j) = round(melt(i).d(j)/8)+12;
            
            %compile thermal forcing data to plot only one thermal forcing 
            %estimate for all icebergs that use the same profile
            TFcoords = [TFcoords; oceanx oceany]; TF = [TF; (Tavg - Tfpavg)'];
            
            clear minrefs oceantemps oceansals oceandepths oceanx oceany timespan Tavg Tfpavg Savg;
        end
        
        %calculate the median thermal forcing for each profile
        [unique_coords,unique_refs,inds] = unique(TFcoords,'rows');
        for j = 1:max(unique_refs)
            %median of thermal forcing
            TFmedian(j) = nanmedian(TF(inds==j));
            
            %create the colormap to show thermal forcing on the map
            if ~isnan(TFmedian(j))
                tempref = round(TFmedian(j)*(2*cmap_mult));
                if tempref < 1
                    temp_map(j,:) = [0 0 0];
                else
                    temp_map(j,:) = Tforcing_cmap(tempref,:);
                end
                clear tempref;
            else
                temp_map(j,:) = [1 1 1];
            end
        end
        %plot the median thermal forcing from each profile on the overview map
        figure(Tm_mapplot);
        scatter(unique_coords(sum(temp_map,2)~=3,1),unique_coords(sum(temp_map,2)~=3,2),16,temp_map((sum(temp_map,2)~=3),:),'filled','o'); hold on;
        
        %plot the median thermal forcing from each profile on the site map
        figure(Tmsitemap);
        scatter(unique_coords(sum(temp_map,2)~=3,1),unique_coords(sum(temp_map,2)~=3,2),16,temp_map((sum(temp_map,2)~=3),:),'filled','o'); hold on;
        sitex = [sitex; unique_coords(sum(temp_map,2)~=3,1)]; sitey = [sitey; unique_coords(sum(temp_map,2)~=3,2)]; 
        clear temp_map;
        
        
        %create an individual plot of temp vs meltrate for the study site
        figure; set(gcf,'position',[850 650 800 400]);
        scatter(melt(i).oceanTavg-melt(i).oceanTfp,100*melt(i).m,2*draft_size,region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
        set(gca,'fontsize',16); grid on;
        xlabel(['Thermal forcing (',char(176),'C above freezing)'],'fontsize',16); ylabel('Melt rate (cm/d)','fontsize',16);
        title(melt(i).dispname); drawnow;
        [f,gof] = fit(melt(i).oceanTavg-melt(i).oceanTfp,melt(i).m,'poly1'); %disp(['Trendline r^2 = ',num2str(gof.rsquare)]);
        if gof.rsquare > 0.5
            saveas(gcf,[figure_path,char(melt(i).name),'-iceberg-oceandata-scatterplot.eps'],'epsc'); saveas(gcf,[figure_path,char(melt(i).name),'-iceberg-oceandata-scatterplot.png'],'png');
        end
        
        %add temp vs meltrate data to composite scatterplot
        figure(Tm_scatterplot);
        if ~isempty(strmatch('Edgeworth',char(leg_names(i)))) %only create a handle for one study site from each region (EAP,EAIS,WAIS,WAP)
            sp(1) = scatter(melt(i).oceanTavg-melt(i).oceanTfp,365*melt(i).m,draft_size,region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
        elseif ~isempty(strmatch('Mertz',char(leg_names(i))))
            sp(2) = scatter(melt(i).oceanTavg-melt(i).oceanTfp,365*melt(i).m,draft_size,region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
        elseif ~isempty(strmatch('Thwaites',char(leg_names(i))))
            sp(3) = scatter(melt(i).oceanTavg-melt(i).oceanTfp,365*melt(i).m,draft_size,region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
        elseif ~isempty(strmatch('Cadman',char(leg_names(i))))
            sp(4) = scatter(melt(i).oceanTavg-melt(i).oceanTfp,365*melt(i).m,draft_size,region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
        else
            scatter(melt(i).oceanTavg-melt(i).oceanTfp,365*melt(i).m,draft_size,region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
        end
        %         for j = 1:length(melt(i).to)
        %             plot(melt(i).oceanTavg(j)-melt(i).oceanTfp(j,1),365*melt(i).m(j),[marker(i),'k'],'markerfacecolor',depth_cmap(round(melt(i).d(j)),:),'markersize',12); hold on;
        %         end
        
        clear draft_size;
    end
    
    %add iceberg locations to the overview map
    figure(Tm_mapplot); clear draft_map;
    for k = 1:length(melt(i).x)
        draft_map(k,:) = depth_cmap(round(nanmean(melt(i).d(k))),:);
    end
    max_ind = find(melt(i).d == max(melt(i).d)); %only plot a symbol for the deepest-drafted iceberg
    scatter(melt(i).x(max_ind),melt(i).y(max_ind),48,draft_map(max_ind,:),'filled','s','markeredgecolor','k','linewidth',0.5); hold on;
    
    %add iceberg locations to the site map
    figure(Tmsitemap);
    scatter(melt(i).x,melt(i).y,48,draft_map,'filled','s','markeredgecolor','k','linewidth',0.5); hold on;
    clear draft_map;
    
    %format the site map & save
    figure(Tmsitemap);
    sitex = [sitex; melt(i).x]; sitey = [sitey; melt(i).y];
    if sqrt((max(sitex)-min(sitex)).^2 + (max(sitey)-min(sitey)).^2)+10000 < 50000
        set(gca,'xlim',[min(sitex)-5000 max(sitex)+5000],'ylim',[min(sitey)-5000 max(sitey)+5000],...
            'xtick',[(ceil(min(sitex)/1000)*1000-5000):5000:(floor(max(sitex)/1000)*1000+5000)],...
            'xticklabel',[(ceil(min(sitex)/1000)-5):5:(floor(max(sitex)/1000)+5)],...
            'ytick',[(ceil(min(sitey)/1000)*1000-5000):5000:(floor(max(sitey)/1000)*1000+5000)],...
            'yticklabel',[(ceil(min(sitey)/1000)-5):5:(floor(max(sitey)/1000)+5)],...
            'fontsize',16);
    else
        set(gca,'xlim',[min(sitex)-5000 max(sitex)+5000],'ylim',[min(sitey)-5000 max(sitey)+5000],...
            'xtick',[(ceil(min(sitex)/1000)*1000-5000):10000:(floor(max(sitex)/1000)*1000+5000)],...
            'xticklabel',[(ceil(min(sitex)/1000)-5):10:(floor(max(sitex)/1000)+5)],...
            'ytick',[(ceil(min(sitey)/1000)*1000-5000):5000:(floor(max(sitey)/1000)*1000+5000)],...
            'yticklabel',[(ceil(min(sitey)/1000)-5):5:(floor(max(sitey)/1000)+5)],...
            'fontsize',16);
    end
    xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
    %resize vertical dimension to maximize figure window usage
    xlims = get(gca,'xlim'); ylims = get(gca,'ylim'); figpos = get(gcf,'position');
    set(gcf,'position',[figpos(1) figpos(2) figpos(3) (max(ylims)-min(ylims))/(max(xlims)-min(xlims))*figpos(3)]);
    clear im;
    saveas(Tmsitemap,[figure_path,char(region(i)),'_iceberg-oceandata-map.eps'],'epsc'); saveas(Tmsitemap,[figure_path,char(region(i)),'_iceberg-oceandata-map.png'],'png');
    clear xlims ylims; close(Tmsitemap);
end

%format the overview map & save
figure(Tm_mapplot);
%add labels
for i = 1:length(melt)
    if strcmp(char(plot_letters(i)),'f)')
        text(nanmean(melt(i).x)-200000,nanmean(melt(i).y)-100000,char(plot_letters(i)),'fontsize',16);
    elseif strcmp(char(plot_letters(i)),'a)')
        text(nanmean(melt(i).x)-200000,nanmean(melt(i).y)+100000,char(plot_letters(i)),'fontsize',16);
    elseif strcmp(char(plot_letters(i)),'b)') || strcmp(char(plot_letters(i)),'c)') || strcmp(char(plot_letters(i)),'e)')
        text(nanmean(melt(i).x)-200000,nanmean(melt(i).y)-100000,char(plot_letters(i)),'fontsize',16);
    else
        text(nanmean(melt(i).x)+100000,nanmean(melt(i).y),char(plot_letters(i)),'fontsize',16);
    end
end
%label axes of map
set(gca,'xlim',[-28e5 28e5],'xtick',[-32e5:8e5:32e5],'xticklabel',[-3200:800:3200],...
    'ylim',[-24e5 24e5],'ytick',[-24e5:8e5:24e5],'yticklabel',[-2400:800:2400],'fontsize',16); grid off;
xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
%add polar stereo coordinates
graticuleps(-50:-5:-90,-180:30:180);
text(0,6.5e5,'85^oS','fontsize',16); text(0,12.0e5,'80^oS','fontsize',16); text(0,17.5e5,'75^oS','fontsize',16); text(0,23.0e5,'70^oS','fontsize',16);
text(-16.5e5,25.25e5,'-30^oE','fontsize',16); text(12.5e5,25.25e5,'30^oE','fontsize',16); 
%add color legend for depth
rectangle('position',[-26.5e5 -23.75e5 16e5 10e5],'facecolor','w','edgecolor','k');
rectangle('position',[-26.05e5 -22.8e5 2.1e5 length(depth_cmap)*1000+0.1e5],'facecolor','k','edgecolor','k');
for k = 1:length(depth_cmap)
    plot([-26e5 -24e5],[-22.75e5+k*1000 -22.75e5+k*1000],'-','color',depth_cmap(end-(k-1),:)); hold on;
end
text(-23.25e5,-22.75e5+k*1000,'0 m','fontsize',16); 
text(-23.25e5,-22.75e5+(k-200)*1000,'200 m','fontsize',16); 
text(-23.25e5,-22.75e5+(k-800)*1000,'800 m','fontsize',16);
%add color legend for thermal forcing
rectangle('position',[-17.30e5 -22.8e5 2.1e5 length(Tforcing_cmap)*1000+0.1e5],'facecolor','k','edgecolor','k');
for k = 1:length(Tforcing_cmap)
    plot([-17.25e5 -15.25e5],[-22.75e5+k*1000 -22.75e5+k*1000],'-','color',Tforcing_cmap(k,:)); hold on;
end
text(-14.5e5,-22.75e5+k*1000,[num2str((length(Tforcing_cmap)/(2*cmap_mult))),char(176),'C'],'fontsize',16); 
text(-14.5e5,-22.75e5+(length(Tforcing_cmap)/2)*1000,[num2str((length(Tforcing_cmap)/(2*cmap_mult))/2),char(176),'C'],'fontsize',16); 
text(-14.5e5,-22.75e5,[num2str(0),char(176),'C'],'fontsize',16);
colormap(gca,im_cmap);%make sure the image colormap didn't get accidentally altered
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
saveas(Tm_mapplot,[figure_path,'Antarctic-iceberg-oceandata-map.eps'],'epsc'); saveas(Tm_mapplot,[figure_path,'Antarctic-iceberg-oceandata-map.png'],'png');


%label the scatterplot & save
figure(Tm_scatterplot);  
set(gca,'fontsize',16); grid on;
%uncomment next 4 lines if you use the plot function to specify symbol colors as a function of draft
% for k = 1:length(depth_cmap)
%     plot([0.1 0.25],[47.5-(k/100) 47.5-(k/100)],'-','color',depth_cmap(k,:)); hold on;
% end
% text(0.275,47,'0 m','fontsize',16); text(0.275,45,'200 m','fontsize',16); text(0.275,39.5,'750 m','fontsize',16);
%next 9 lines should be used if scatterplot function specifies symbol colors as a function of region & symbol size as a function of draft
ylims = get(gca,'ylim');
rectangle('position',[0.425 max(ylims) - 0.06*((depth_cutoff-50)/150*1.15)*(range(ylims)) 0.3 0.06*((depth_cutoff-50)/150*1.05)*(range(ylims))],'facecolor','w','edgecolor','k');
for j = 1:1:(depth_cutoff-50)/150
    draft_size(j) = round((50+((j-1)*150))/8)+12;
    yloc(j) = max(ylims) - 0.06*j*(range(ylims)) - 0.01*(range(ylims));
    text(0.525,yloc(j),[num2str((50+((j-1)*150))),' m'],'fontsize',16);
end
scatter(repmat(0.475,size(yloc)),yloc,draft_size,'w','filled','s','markeredgecolor','k'); hold on;
sp_leg = legend(sp,'EAP','EAIS','WAIS','WAP'); set(sp_leg,'location','northwest');
xlabel(['Thermal forcing (',char(176),'C above freezing)'],'fontsize',16); ylabel('Melt rate (m/yr)','fontsize',16);
saveas(Tm_scatterplot,[figure_path,'Antarctic-iceberg-meltrate-temp-depth-scatterplots.eps'],'epsc'); saveas(Tm_scatterplot,[figure_path,'Antarctic-iceberg-meltrate-temp-depth-scatterplots.png'],'png');
clear ylims yloc draft_size;

%resave the data structure
save([iceberg_path,'Antarctic-icebergmelt-comparison.mat'],'melt','-v7.3');


% %make maps showing the average temperature in the upper 200m & the maximum
% %temperature
% Temp_mapplot = figure; set(gcf,'position',[450 50 800 800]); Temp_cmap = colormap(hot(100));
% imagesc(IM.x,IM.y,IM.z); colormap gray; hold on; axis xy equal; colormap(gray(10001));
% cmap_ref = [1 1 2 2 2 2 2 3 3 4 4 4 4 4 4];
% for i = 1:length(melt)
%     if ~isempty(melt(i).oceant)
%         if size(melt(i).oceanx,1) > 1
%             for j = 1:size(melt(i).oceanx,1)
%                 temp(j,1) = nanmean(melt(i).oceanT(melt(i).oceand(:,j)<200,j));
%                 freeze_temp(j,1) = nanmean(((-5.73*10^-2).*melt(i).oceanS(melt(i).oceand(:,j)<200,j)) + (8.32*10^-2) - ((7.61*10^-4).*melt(i).oceand(melt(i).oceand(:,j)<200,j)));
%                 if ~isnan(freeze_temp(j,1)) && ~isnan(temp(j,1))
%                     plot(melt(i).oceanx(j),melt(i).oceany(j),'ok','markerfacecolor',Temp_cmap(round((temp(j,1)-freeze_temp(j,1)+0.5)*30),:)); hold on;
%                 end
%             end
%         else
%             for j = 1:size(melt(i).oceanx,2)
%                 temp(j,1) = nanmean(melt(i).oceanT(melt(i).oceand(:,j)<200,j));
%                 freeze_temp(j,1) = nanmean(((-5.73*10^-2).*melt(i).oceanS(melt(i).oceand(:,j)<200,j)) + (8.32*10^-2) - ((7.61*10^-4).*melt(i).oceand(melt(i).oceand(:,j)<200,j)));
%                 if ~isnan(freeze_temp(j,1)) && ~isnan(temp(j,1))
%                     plot(melt(i).oceanx(j),melt(i).oceany(j),'ok','markerfacecolor',Temp_cmap(round((temp(j,1)-freeze_temp(j,1)+0.5)*25),:)); hold on;
%                 end
%             end
%         end
%         drawnow;
%     end
% end

