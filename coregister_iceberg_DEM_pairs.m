function [DEM1,DEM2] = coregister_iceberg_DEM_pairs(varargin)
% Function to vertically co-register iceberg DEMS
% Ellyn Enderlin (ellynenderlin@boisestate.edu)
% Last edited: 09 Nov. 2022
% 
% INPUTS:   DEM1            structure variable containing earlier DEM info
%           DEM2            structure variable containing later DEM info
%           IM1             structure variable containing earlier
%                               orthoimage info
%           IM2             structure variable containing later orthoimage
%                               info
%           geography       binary specification of polar region (0 = Greenland, 1 = Antarctic)
%           region_abbrev   region abbreviation used in image files
%           dir_TMD         directory path to the TMD folder ('*/TMD/')
%           dir_output      directory where all output files will be placed
%
% OUTPUTS:  DEM1            DEM1 structure variable with new fields
%           DEM2            DEM2 structure variable with new fields
%
% Calls the following external functions:
%   - ps2wgs.m
%   - tmd_tide_pred.m
%   - extract_geoid_elevs.m
%   - calculate_DEM_offset.m

%assign variable names & specify options based on inputs
DEM1 = varargin{1}; DEM2 = varargin{2};
IM1 = varargin{3}; IM2 = varargin{4}; 
geography = varargin{5}; region_abbrev = varargin{6};
dir_output = varargin{7};
if nargin == 8
    dir_TMD = varargin{8};
    tidemodel_flag = 1;
else
    tidemodel_flag = 0;
end

% ----------plot the tidal model map & select the nearest neighboring grid cell----------
disp('Creating a text file specifying the area over which geoid elevations (& optional tidal heights) should be pulled...');
%convert to lat-lon (mapxy inputs in km, not m like ps2wgs)
x_min = min([DEM1.x DEM2.x]); x_max = max([DEM1.x DEM2.x]);
y_min = min([DEM1.y DEM2.y]); y_max = max([DEM1.y DEM2.y]);
if geography == 0
    [ROIlon, ROIlat] = ps2wgs([x_min x_max x_max x_min],[y_min y_min y_max y_max],'StandardParallel',70,'StandardMeridian',-45); %Greenland PS
elseif geography == 1
    [ROIlon, ROIlat] = ps2wgs([x_min x_max x_max x_min],[y_min y_min y_max y_max],'StandardParallel',-71,'StandardMeridian',0); %Antarctic PS
end
C=[min(ROIlat) min(ROIlon) max(ROIlat) max(ROIlon)]; %ULlat ULlon LRlat LRlon
dlmwrite([dir_output,region_abbrev,'_corner_WGScoords.txt'],C,'\t');

% Create tidal coordinates file if it does not already exist
if tidemodel_flag == 1 && ~isfile([dir_output,region_abbrev,'_tidal_coords.txt'])
    %create tidal map around the region of interest
    early_time = datenum(DEM1.time,'yyyymmddHHMMSS');
    late_time = datenum(DEM2.time,'yyyymmddHHMMSS');
    SD_time = (late_time+early_time)/2;
    lonvec=[min(ROIlon)-0.1:0.05:max(ROIlon)+0.1]; latvec=[min(ROIlat)-0.1:0.05:max(ROIlat)+0.1]';
    lon=repmat(lonvec,length(latvec),1); lat=repmat(latvec,1,length(lonvec));
    disp('Extracting regional tidal height snapshot map');
    cd([dir_TMD,'DATA/']);
    [TS,~] = tmd_tide_pred('Model_CATS2008',SD_time,lat,lon,'z'); %regional snapshot of tidal heights
    
    % Plot results
    figure(2); clf;
    hold on; set(gcf,'position',[50 50 700 700]);
    imagesc(lonvec,latvec,TS);
    axis xy equal; colormap jet; colorbar; set(gca,'fontsize',14,'linewidth',2);
    xlabel('lon'); ylabel('lat');
    %overlay DEM outlines
    plot([ROIlon;ROIlon(1)],[ROIlat;ROIlat(1)],'-w','linewidth',3);
    plot(mean(ROIlon,'omitnan'),mean(ROIlat,'omitnan'),'xw','linewidth',3,'markersize',20);
    [tide,~] = tmd_tide_pred('Model_CATS2008',SD_time,mean(ROIlat,'omitnan'),mean(ROIlon,'omitnan'),'z');
    disp(['Tidal height snapshot value = ',num2str(tide),'m']);
    prompt = 'Specify coordinates for tidal height estimates. Do the coordinates need to differ from the centroid (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        disp('To change coordinates, identify new coordinates using tick marks and type center_lon=XX; center_lat=YY; dbcont ');
        keyboard
    else
        center_lon=mean(ROIlon,'omitnan'); center_lat=mean(ROIlat,'omitnan');
    end
    close all; drawnow;
    
    dlmwrite([dir_output,region_abbrev,'_tidal_coords.txt'],[center_lon center_lat],'\t');
end

% ----------load coordinates and pull geoid elevations----------
disp('Calculating geoid elevations...');
cd(dir_output);
C = dlmread([dir_output,region_abbrev,'_corner_WGScoords.txt'],'\t');
ULlat = C(:,1); ULlon = C(:,2); LRlat = C(:,3); LRlon = C(:,4);

% pull geoid elevations over a ~500 m-resolution grid
[geoid] = extract_geoid_elevs(ULlat,ULlon,LRlat,LRlon,geography,region_abbrev,dir_output); 
% create a mesh grid from coordinate vectors
[DEM1.x_grid,DEM1.y_grid] = meshgrid(DEM1.x,DEM1.y);
[DEM2.x_grid,DEM2.y_grid] = meshgrid(DEM2.x,DEM2.y);
% interpolate geoid elevations to DEM pixel locations
DEM1.geoid_z = griddata(geoid.x,geoid.y,geoid.z,DEM1.x_grid,DEM1.y_grid); 
DEM2.geoid_z = griddata(geoid.x,geoid.y,geoid.z,DEM2.x_grid,DEM2.y_grid); 
% calculate difference between ellipsoidal and geoid heights
DEM1.z_elpsd_adjust = DEM1.z - DEM1.geoid_z; 
DEM2.z_elpsd_adjust = DEM2.z - DEM2.geoid_z; 
% resave DEM matfiles with geoid variables
DEM = DEM1; save([dir_output,DEM1.filename,'.mat'],'DEM','-v7.3');
clear DEM; DEM = DEM2; save([dir_output,DEM2.filename,'.mat'],'DEM','-v7.3');
disp('Resaved DEMs w/ orthometric heights');
clear DEM geoid;

% ----------extract co-registration information----------
disp('Coregistering DEMs...');
% extract center tidal data coordinates
center_coords = [nanmean([ULlon LRlon]), nanmean([ULlat LRlat])];
center_lon = center_coords(1); center_lat = center_coords(2); 

% calculate offset between DEMs
if tidemodel_flag == 0
    [DEM1,DEM2] = calculate_DEM_offset(DEM1,DEM2,IM1,IM2,dir_output,center_lon,center_lat);
else
    [DEM1,DEM2] = calculate_DEM_offset(DEM1,DEM2,IM1,IM2,dir_output,center_lon,center_lat,dir_TMD);
end
close all;

% plot DEMs
DEM1.z_elpsd_adjust(DEM1.z_elpsd_adjust<0) = 0; DEM2.z_elpsd_adjust(DEM2.z_elpsd_adjust<0) = 0; %remove elevations < 0 for plotting purposes
figure1 = figure; set(figure1,'position',[0 600 700 400]);
imagesc(DEM1.x,DEM1.y,DEM1.z_elpsd_adjust); axis xy equal; set(gca,'clim',[min(DEM1.z_elpsd_adjust(~isnan(DEM1.z_elpsd_adjust))) min(DEM1.z_elpsd_adjust(~isnan(DEM1.z_elpsd_adjust)))+40]); colormap(gca,'jet'); colorbar; grid on; hold on;
set(gca,'xtick',[min(DEM1.x):500:max(DEM1.x)],'xticklabel',[min(DEM1.x)/1000:0.5:max(DEM1.x)/1000],...
    'ytick',[min(DEM1.y):500:max(DEM1.y)],'yticklabel',[min(DEM1.y)/1000:0.5:max(DEM1.y)/1000]);
figure2 = figure; set(figure2,'position',[675 600 700 400]);
imagesc(DEM2.x,DEM2.y,DEM2.z_elpsd_adjust); axis xy equal; set(gca,'clim',[min(DEM2.z_elpsd_adjust(~isnan(DEM2.z_elpsd_adjust))) min(DEM2.z_elpsd_adjust(~isnan(DEM2.z_elpsd_adjust)))+40]); colormap(gca,'jet'); colorbar; grid on; hold on;
set(gca,'xtick',[min(DEM1.x):500:max(DEM1.x)],'xticklabel',[min(DEM1.x)/1000:0.5:max(DEM1.x)/1000],...
    'ytick',[min(DEM1.y):500:max(DEM1.y)],'yticklabel',[min(DEM1.y)/1000:0.5:max(DEM1.y)/1000]);

% plot images
figure3 = figure; set(figure3,'position',[0 100 700 400]);
imagesc(IM1.x,IM1.y,IM1.z); axis xy equal; colormap(gca,'gray'); grid on; set(gca,'clim',[0 200]); hold on;
set(gca,'xtick',[min(DEM1.x):500:max(DEM1.x)],'xticklabel',[min(DEM1.x)/1000:0.5:max(DEM1.x)/1000],...
    'ytick',[min(DEM1.y):500:max(DEM1.y)],'yticklabel',[min(DEM1.y)/1000:0.5:max(DEM1.y)/1000]);
figure4 = figure; set(figure4,'position',[675 100 700 400]);
L = imagesc(IM2.x,IM2.y,IM2.z); axis xy equal; colormap(gca,'gray'); grid on; set(gca,'clim',[0 200]); hold on;
set(gca,'xtick',[min(DEM1.x):500:max(DEM1.x)],'xticklabel',[min(DEM1.x)/1000:0.5:max(DEM1.x)/1000],...
    'ytick',[min(DEM1.y):500:max(DEM1.y)],'yticklabel',[min(DEM1.y)/1000:0.5:max(DEM1.y)/1000]);

disp('Advance to the next step');
disp('------------------------');

end

