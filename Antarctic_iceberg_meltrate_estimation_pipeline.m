%% Antarctic Iceberg Meltrate Estimation Pipeline
% Ellyn Enderlin & Rainey Aberle, Fall 2021
% 
% This file serves as a wrapper for the iceberg melt estimation codes.
% Run the following codes in order to estimate iceberg melt rates 
% from digital elevation models.
%
% Files required for each DEM:
%   - "[region_abbrev]_[DEM_time]-DEM.tif"
%   - "[region_abbrev]_[DEM_time]-orthoimage.tif"
%
% To do: 
%   - Edit function descriptions
%   - Make sure DEMs and IMs are saved at the end of every section (and
%       only then)
%   - Ensure text files are exported with headers

%% specify site- and directory-specific information

clearvars; close all;

% ----------MODIFY VARIABLES BELOW----------
% site name and region abbreviation used for file names
region_name = 'Edgeworth-LarsenA'; 
region_abbrev = 'LA';
% path to Iceberg-melt-code in directory
dir_code = '/Users/icebergs/iceberg-melt/'; 
% path to TMD folder in directory
dir_TMD = [dir_code,'general/TMD/'];
% path to DEM files in directory
dir_DEM = [dir_code,region_name,'/DEMs/'];
% path to outputs folder (where you would like all outputs saved)
% dir_output = [dir_code,'../meltrate_outputs/'];
dir_output = dir_DEM;
% DEM time stamps (DEM1 = earlier, DEM2 = later) used in file names
% (YYYYMMDDhhmmss)
DEM1.time = '20191012101214'; DEM2.time = '20200927092713'; 

% ----------INITIAL SET-UP----------
% General DEM filenames (no suffixes) - Ensure filenames match this format
DEM1.filename = [region_abbrev,'_',DEM1.time];
DEM2.filename = [region_abbrev,'_',DEM2.time];

% add paths to necessary functions and datasets
addpath(dir_code);
addpath([dir_code,'functions/']);
addpath([dir_code,'general/']);
addpath([dir_code,'general/TMD/']);
addpath([dir_code,'general/TMD/FUNCTIONS/']);
addpath([dir_code,'general/RACMO2.3_Antarctica/']);

%% 1. Convert DEMs & accompanying images to mat-files
cd_to_sitedir = ['cd ',dir_output]; eval(cd_to_sitedir);

% create matfiles from geotiffs for the earlier date DEM
[DEM1,IM1] = convert_ASP_tifs_to_matfiles(DEM1,dir_DEM,dir_output);
% Plot
figure(1); set(gcf,'position',[0 200 800 800]);
ax1 = subplot(2,2,1);  
    imagesc(ax1,IM1.x/10^3,IM1.y/10^3,IM1.z); axis xy equal; colormap gray;
    ax1.YLabel.String = 'Northing (km)'; ax1.FontSize = 14; title('DEM1');
ax3 = subplot(2,2,3);
    imagesc(ax3,DEM1.x/10^3,DEM1.y/10^3,DEM1.z); axis xy equal; colormap hot;
    ax3.YLabel.String = 'Northing (km)'; ax3.XLabel.String='Easting (km)'; ax3.FontSize = 14;

% create matfiles from geotiffs for the earlier date DEM
[DEM2,IM2] = convert_ASP_tifs_to_matfiles(DEM2,dir_DEM,dir_output);
% Plot
ax2 = subplot(2,2,2); 
    imagesc(ax2,IM2.x./10^3,IM2.y./10^3,IM2.z); axis xy equal; colormap gray;
    ax2.YLabel.String = 'Northing (km)'; ax2.FontSize = 14; title('DEM2'); 
ax4 = subplot(2,2,4); 
    imagesc(ax4, DEM2.x./10^3,DEM2.y./10^3,DEM2.z); axis xy equal; colormap hot;
    ax4.YLabel.String = 'Northing (km)'; ax4.XLabel.String='Easting (km)'; ax4.FontSize = 14;
    %colorbar;

clear ax*;

disp('Advance to the next step');
disp('------------------------');

%% 2. Vertically coregister DEMs to account for satellite uncertainty & tidal change between dates
cd_to_sitedir = ['cd ',dir_output]; eval(cd_to_sitedir);
close all;

% Option to reload MAT files from previous step if previously completed but
% variables are not in workspace
reload = input('Would you like to load DEM and IM variables from file (y/n)?','s');
if strcmp(reload,'y')
    DEM1 = load([dir_output,region_abbrev,'_',DEM1.time,'-DEM.mat']).DEM;
    IM1 = load([dir_output,region_abbrev,'_',DEM1.time,'-orthoimage.mat']).IM;
    DEM2 = load([dir_output,region_abbrev,'_',DEM2.time,'-DEM.mat']).DEM;
    IM2 = load([dir_output,region_abbrev,'_',DEM2.time,'-orthoimage.mat']).IM;
end
disp('Files loaded');

[DEM1,DEM2] = coregister_iceberg_DEM_pairs(DEM1,DEM2,IM1,IM2,region_abbrev,dir_TMD,dir_output);

%% 3. Identify and save iceberg coordinates
cd_to_sitedir = ['cd ',dir_output]; eval(cd_to_sitedir);
close all; 

%create a date-specific directory as necessary
if not(isfolder([dir_output,'/',DEM1.time,'-',DEM2.time,'/']))
    mkdir([dir_output,'/',DEM1.time,'-',DEM2.time,'/']) % make models folder if it doesn't exist
    disp('date-specific folder created.');
end

% Option to reload MAT files from previous step if previously completed but
% variables are not in workspace
reload = input('Would you like to load DEM and IM variables from file (y/n)?','s');
if strcmp(reload,'y')
    DEM1 = load([dir_output,region_abbrev,'_',DEM1.time,'-DEM.mat']).DEM;
    IM1 = load([dir_output,region_abbrev,'_',DEM1.time,'-orthoimage.mat']).IM;
    DEM2 = load([dir_output,region_abbrev,'_',DEM2.time,'-DEM.mat']).DEM;
    IM2 = load([dir_output,region_abbrev,'_',DEM2.time,'-orthoimage.mat']).IM;
    disp('files loaded');
end

[PSx_early,PSy_early,PSx_late,PSy_late] = track_icebergs(DEM1,DEM2,IM1,IM2,dir_output); 

%% 4. Extract elevation change, convert to volume change, then estimate melt rate
cd_to_sitedir = ['cd ',dir_output]; eval(cd_to_sitedir);
close all; 

% Select step to run in estimate_iceberg_meltrates function:
%   1 = Estimate elevation change for each iceberg
%   2 = UPDATE individual icebergs and/or remove icebergs with anomalous melt rate estimates
step_no = 1; 

% Option to reload MAT files from previous step if previously completed but
% variables are not in workspace
reload = input('Would you like to load DEM and IM variables from file (y/n)?','s');
if strcmp(reload,'y')
    DEM1 = load([dir_output,region_abbrev,'_',DEM1.time,'-DEM.mat']).DEM;
    IM1 = load([dir_output,region_abbrev,'_',DEM1.time,'-orthoimage.mat']).IM;
    DEM2 = load([dir_output,region_abbrev,'_',DEM2.time,'-DEM.mat']).DEM;
    IM2 = load([dir_output,region_abbrev,'_',DEM2.time,'-orthoimage.mat']).IM;
    disp('DEM and IM variables loaded');
end

% estimate iceberg melt rates
[DEM1,DEM2] = estimate_iceberg_meltrates(DEM1,DEM2,IM1,IM2,dir_output,dir_code,region_abbrev,region_name,step_no);
