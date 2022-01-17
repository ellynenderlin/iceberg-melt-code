function [geoid] = extract_Antarctic_geoid_elevs(ULlat,ULlon,LRlat,LRlon,region_abbrev,dir_output)
% Function to extract Antarctic geoid elevations over a given region
% Ellyn Enderlin & Mariama Dryak 
% Slightly reformatted by Rainey Aberle (Fall 2021)
% 
% INPUTS:   ULlat           latitude of image upper left corner 
%           ULlon           longitude of image upper left corner
%           LRlat           latitude of image lower right corner
%           LRlon           longitude of image lower right corner
%           region_abbrev   region abbreviation used in image files
%           dir_output      directory where all output files will be placed
%
% OUTPUTS:  DEM1            DEM1 structure variable with new fields
%           DEM2            DEM2 structure variable with new fields
%
% Calls the following external functions:
%   - wgs2ps.m

%create a mat-file containing geoid elevations over the ROI containing your DEMs on a 500 m-resolution grid
%cd_root_dir = ['cd ',root_dir,region_name]; eval(cd_root_dir);
SP = -71; SM = 0; %Antarctic Polar Stereo standard parallel & meridian
if ULlon < 0
    ULlon = 360+ULlon; LRlon = 360+LRlon;
end
max_WGSlat = max([ULlat;LRlat])+1; min_WGSlat = min([ULlat;LRlat])-1;
max_WGSlon = max([ULlon;LRlon])+1; min_WGSlon = min([ULlon;LRlon])-1;
WGSlat = min_WGSlat:0.01:max_WGSlat;
WGSlon = min_WGSlon:0.01:max_WGSlon;
[ylon,ylat] = meshgrid(WGSlon,WGSlat);
for i = 1:size(ylat,1)
    for j = 1:size(ylat,2)
        geoid_z(i,j) = geoidheight(ylat(i,j),ylon(i,j));
    end
end
ylon = ylon-360;
[psx_grid,psy_grid] = wgs2ps(ylon,ylat,'StandardParallel',SP,'StandardMeridian',SM);

geoid.x = psx_grid; geoid.y = psy_grid; geoid.z = geoid_z;
geoid.lon = WGSlon-360; geoid.lat = WGSlat;
save([dir_output,region_abbrev,'_geoid_heights.mat'],'geoid','-v7.3');
disp('geoid heights saved');

end