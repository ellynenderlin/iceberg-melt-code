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
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code','/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean');
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/AntarcticMappingTools');

%specify paths & file names for data
CTD_path = '/Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/CTD_Antarctica/';
CTD_data = [CTD_path,'Antarctic-ocean-data.mat'];
RACMO_path = '/Users/ellynenderlin/Research/miscellaneous/RACMO2.3_Antarctica/';
iceberg_path = '/Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt/';
figure_path = [iceberg_path,'figures/'];

%specify study site names
region = [{'Ronne'},{'Filchner'},{'Amery'},{'Totten'},{'Mertz'},{'Edgeworth-LarsenA'},{'Crane-LarsenB'},...
    {'Thwaites'},{'Ferrigno-Eltanin'},{'Seller-Bugge'},{'Heim-Marguerite'},{'Widdowson-Biscoe'},{'Cadman-Palmer'},{'Blanchard-Danco'},{'Leonardo-Danco'}];
leg_names = [{'Edgeworth'},{'Crane'},{'Ronne'},{'Filchner'},{'Polar Times'},{'Totten'},{'Mertz'},{'Thwaites'},{'Ferrigno'},{'Seller'},{'Heim'},{'Widdowson'},{'Cadman'},{'Blanchard'},{'Leonardo'}];
leg_ref = [9,10,11,12,13,14,15,8,7,6,5,4,3,2,1]; %arrange legend in alphabetical order

%specify plot params
marker = ['s','s','s','s','s','s','s','s','s','s','s','s','s','s','s']; %modify if you want to change symbols to indicate something about data (E vs W for example)
plot_letters = [{'k)'},{'l)'},{'m)'},{'n)'},{'o)'},{'i)'},{'j)'},{'h)'},{'g)'},{'f)'},{'e)'},{'d)'},{'c)'},{'b)'},{'a)'}]; %plot letters for sites to be used in geographically-arranged subplots
plot_loc = [6,8,10,12,14,2,4,15,13,11,9,7,5,3,1];
region_colors = [184,225,134; 184,225,134; 184,225,134; 184,225,134; 184,225,134; 39,100,25; 39,100,25;...
    241,182,218; 241,182,218; 197,27,125; 197,27,125; 197,27,125; 197,27,125; 197,27,125; 197,27,125]./255; 
cmap_add = 2.5; cmap_mult = 100;
Temp_cmap = cmocean('thermal',cmap_add*200); 
highmelt_cmap = cmocean('amp',100);
Tforcing_cmap = cmocean('amp',600);

%specify generic variables
rho_sw = 1026; %sea water density in kg m^-3
depth_cutoff = 800; %maximum depth of icebergs & therefore ocean observation data of interest in meters
years = [2011.75 2022.25]; year_ticks = [2013:2:2022]; %approximate date range for plots

%specify the buffer region to search for ocean data around the icebergs
buffer = 100000; %m

%load Antarctic image to plot as background for a map
Antarctic_map = 'LIMA';
cd(['/Users/ellynenderlin/Research/miscellaneous/',Antarctic_map]);
if strcmp('LIMA',Antarctic_map)
    [A,S] = readgeoraster('00000-20080319-092059124.tif');
    IM.x = S.XWorldLimits(1):S.SampleSpacingInWorldX:S.XWorldLimits(2);
    IM.y = S.YWorldLimits(2):-S.SampleSpacingInWorldY:S.YWorldLimits(1);
    IM.z=single(A); IM.z = IM.z./255;
else %assume REMA
    IM.x =S.XWorldLimits(1)+0.5*S.CellExtentInWorldX:S.CellExtentInWorldX:S.XWorldLimits(2)-0.5*S.CellExtentInWorldX;
    IM.y = S.YWorldLimits(2)-0.5*S.CellExtentInWorldY:-S.CellExtentInWorldY:S.YWorldLimits(1)+0.5*S.CellExtentInWorldY;
    IM.z=single(A);
end
clear A S;

%navigate to the iceberg directory as the default workspace
cd(iceberg_path);
close all;

%% organize ocean data & pair with remotely-sensed iceberg melt data (run once to create Antarctic-icebergmelt-comparison.mat)
% If you rerun this section, you will overwrite
% Antarctic-icebergmelt-comparison.mat... rerun if you want to change the
% size of the "buffer" search window for ocean data
disp('Compiling iceberg melt & nearby ocean observations for each site');

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

%% interpolate profiles to a uniform, site-specific profile for each site
close all; warning off;
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


%% compute thermal forcing 
close all; drawnow;
disp('Solving for thermal forcing');

%reload compiled data as needed
if ~exist('melt')
    load([iceberg_path,'Antarctic-icebergmelt-comparison.mat']);
end

%set-up a standard depth profile
tempref = [];
for i = 1:size(region,2)
    disp_names(i) = {strjoin([cellstr(plot_letters(i)),' ',cellstr(leg_names(i))])};
    
    %identify the maximum iceberg draft
    max_ind = find(melt(i).d == max(melt(i).d)); %identify the deepest iceberg
    
    %find appropriate ocean data & extract thermal forcing estimates
    if ~isempty(melt(i).oceant)
        
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
            clear minrefs oceantemps oceansals oceandepths oceanx oceany timespan Tavg Tfpavg Savg;
        end
        
    end
    
end

%format the overview map & save
figure(Tm_mapplot);
%add labels
for i = 1:length(melt)
    if strcmp(char(plot_letters(i)),'f)')
        text(nanmean(melt(i).x)-100000,nanmean(melt(i).y)-100000,char(plot_letters(i)),'fontsize',12);
    elseif strcmp(char(plot_letters(i)),'a)')
        text(nanmean(melt(i).x)-100000,nanmean(melt(i).y)+75000,char(plot_letters(i)),'fontsize',12);
    elseif strcmp(char(plot_letters(i)),'b)') || strcmp(char(plot_letters(i)),'c)') || strcmp(char(plot_letters(i)),'e)')
        text(nanmean(melt(i).x)-100000,nanmean(melt(i).y)-50000,char(plot_letters(i)),'fontsize',12);
    else
        text(nanmean(melt(i).x)+75000,nanmean(melt(i).y),char(plot_letters(i)),'fontsize',12);
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
rectangle('position',[-26.5e5 -23.5e5 28e5 15e5],'facecolor','w','edgecolor','k'); xlims = get(gca,'xlim'); ylims = get(gca,'ylim');

%add color & size legends for iceberg data
%sizes
scatter(min(xlims)+0.4*(max(xlims)-min(xlims)),min(ylims)+0.07*(max(ylims)-min(ylims)),round(100/5+10),'w','filled','s','markeredgecolor','k'); text(min(xlims)+0.425*(max(xlims)-min(xlims)),min(ylims)+0.07*(max(ylims)-min(ylims)),'100 m','fontsize',16); %scaling y-offset = 0.16*(max(ylims)-min(ylims))
scatter(min(xlims)+0.4*(max(xlims)-min(xlims)),min(ylims)+0.15*(max(ylims)-min(ylims)),round(300/5+10),'w','filled','s','markeredgecolor','k'); text(min(xlims)+0.425*(max(xlims)-min(xlims)),min(ylims)+0.15*(max(ylims)-min(ylims)),'300 m','fontsize',16); %scaling y-offset = 0.11*(max(ylims)-min(ylims))
scatter(min(xlims)+0.4*(max(xlims)-min(xlims)),min(ylims)+0.23*(max(ylims)-min(ylims)),round(500/5+10),'w','filled','s','markeredgecolor','k'); text(min(xlims)+0.425*(max(xlims)-min(xlims)),min(ylims)+0.98*0.23*(max(ylims)-min(ylims)),'500 m','fontsize',16); %scaling y-offset = 0.05*(max(ylims)-min(ylims))
text(min(xlims)+0.40*(max(xlims)-min(xlims)),min(ylims)+0.305*(max(ylims)-min(ylims)),'iceberg','fontsize',16,'fontweight','bold');
text(min(xlims)+0.415*(max(xlims)-min(xlims)),min(ylims)+0.275*(max(ylims)-min(ylims)),'draft','fontsize',16,'fontweight','bold');
%colors
for k = 1:length(highmelt_cmap)
    plot([min(xlims)+0.20*(max(xlims)-min(xlims)) min(xlims)+0.25*(max(xlims)-min(xlims))],...
        [min(ylims)+0.245*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap)) min(ylims)+0.245*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap))],...
        '-','linewidth',2*((max(ylims)-min(ylims))/(max(xlims)-min(xlims))),'color',highmelt_cmap((length(highmelt_cmap)+1)-k,:));
end
text(min(xlims)+0.26*(max(xlims)-min(xlims)),min(ylims)+0.245*(max(ylims)-min(ylims)),'50 m yr^{-1}','fontsize',16);
text(min(xlims)+0.26*(max(xlims)-min(xlims)),min(ylims)+0.245*(max(ylims)-min(ylims))-(length(highmelt_cmap)*4/5)*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap)),'10 m yr^{-1}','fontsize',16);
text(min(xlims)+0.26*(max(xlims)-min(xlims)),min(ylims)+0.245*(max(ylims)-min(ylims))-(0.98*length(highmelt_cmap)*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap))),'<1 m yr^{-1}','fontsize',16);
text(min(xlims)+0.23*(max(xlims)-min(xlims)),min(ylims)+0.305*(max(ylims)-min(ylims)),'iceberg','fontsize',16,'fontweight','bold');
text(min(xlims)+0.22*(max(xlims)-min(xlims)),min(ylims)+0.275*(max(ylims)-min(ylims)),'melt rate','fontsize',16,'fontweight','bold');

% %add color legend for depth
% rectangle('position',[-26.05e5 -22.8e5 2.1e5 length(depth_cmap)*1000+0.1e5],'facecolor','k','edgecolor','k');
% for k = 1:length(depth_cmap)
%     plot([-26e5 -24e5],[-22.75e5+k*1000 -22.75e5+k*1000],'-','color',depth_cmap(end-(k-1),:)); hold on;
% end
% text(-23.25e5,-22.75e5+k*1000,'0 m','fontsize',16); 
% text(-23.25e5,-22.75e5+(k-200)*1000,'200 m','fontsize',16); 
% text(-23.25e5,-22.75e5+(k-800)*1000,'800 m','fontsize',16);

% %add color legend for thermal forcing
% rectangle('position',[-17.30e5 -22.8e5 2.1e5 length(Tforcing_cmap)*1000+0.1e5],'facecolor','k','edgecolor','k');
% for k = 1:length(Tforcing_cmap)
%     plot([-17.25e5 -15.25e5],[-22.75e5+k*1000 -22.75e5+k*1000],'-','color',Tforcing_cmap(k,:)); hold on;
% end
% text(-14.5e5,-22.75e5+k*1000,[num2str((length(Tforcing_cmap)/(2*cmap_mult))),char(176),'C'],'fontsize',16); 
% text(-14.5e5,-22.75e5+(length(Tforcing_cmap)/2)*1000,[num2str((length(Tforcing_cmap)/(2*cmap_mult))/2),char(176),'C'],'fontsize',16); 
% text(-14.5e5,-22.75e5,[num2str(0),char(176),'C'],'fontsize',16);
%add color legend for ocean temperature
% rectangle('position',[-17.30e5 -22.8e5 2.1e5 length(Temp_cmap)*1000+0.1e5],'facecolor','k','edgecolor','k');
for k = 1:length(Temp_cmap)
    plot([min(xlims)+0.06*(max(xlims)-min(xlims)) min(xlims)+0.11*(max(xlims)-min(xlims))],...
        [min(ylims)+0.245*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(Temp_cmap)) min(ylims)+0.245*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(Temp_cmap))],...
        '-','color',Temp_cmap((length(Temp_cmap)+1)-k,:)); hold on;
end
text(min(xlims)+0.12*(max(xlims)-min(xlims)),min(ylims)+0.245*(max(ylims)-min(ylims)),[num2str(cmap_add),char(176),'C'],'fontsize',16); 
text(min(xlims)+0.12*(max(xlims)-min(xlims)),min(ylims)+0.245*(max(ylims)-min(ylims))-(length(Temp_cmap)/2)*((0.20*(max(ylims)-min(ylims)))/length(Temp_cmap)),['0',char(176),'C'],'fontsize',16); 
text(min(xlims)+0.12*(max(xlims)-min(xlims)),min(ylims)+0.245*(max(ylims)-min(ylims))-(length(Temp_cmap))*((0.20*(max(ylims)-min(ylims)))/length(Temp_cmap)),['-',num2str(cmap_add),char(176),'C'],'fontsize',16);
text(min(xlims)+0.07*(max(xlims)-min(xlims)),min(ylims)+0.305*(max(ylims)-min(ylims)),'ocean','fontsize',16,'fontweight','bold');
text(min(xlims)+0.075*(max(xlims)-min(xlims)),min(ylims)+0.275*(max(ylims)-min(ylims)),'temp.','fontsize',16,'fontweight','bold');

%save figure
[sorted,inds] = sort(leg_ref); mp_sort = mp(inds);
colormap(gca,im_cmap);%make sure the image colormap didn't get accidentally altered
legmap = legend(mp_sort,[char(disp_names(inds))]); set(legmap,'location','northoutside','fontsize',16,'NumColumns',5); 
legmappos = get(legmap,'position'); set(legmap,'position',[0.05 legmappos(2)+0.05 legmappos(3) legmappos(4)]);
gcapos = get(gca,'position'); set(gca,'position',[gcapos(1) 0.09 gcapos(3) gcapos(4)]);
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





