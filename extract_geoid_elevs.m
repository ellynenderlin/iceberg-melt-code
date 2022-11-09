function [geoid] = extract_geoid_elevs(ULlat,ULlon,LRlat,LRlon,geography,region_abbrev,dir_output)
% Function to extract geoid elevations over a given region
% Ellyn Enderlin (ellynenderlin@boisestate.edu)
% Last edited: 09 Nov. 2022
% 
% INPUTS:   ULlat           latitude of image upper left corner 
%           ULlon           longitude of image upper left corner
%           LRlat           latitude of image lower right corner
%           LRlon           longitude of image lower right corner
%           geography       binary specification of polar region (0 = Greenland, 1 = Antarctic)
%           region_abbrev   region abbreviation used in image files
%           dir_output      directory where all output files will be placed
%
% OUTPUTS:  DEM1            DEM1 structure variable with new fields
%           DEM2            DEM2 structure variable with new fields
%
% Calls the following external functions:
%   - wgs2ps.m

%specify polar projection parameters
if geography == 0
    PSparallel = 70; PSmeridian = -45; %Greenland PS standard parallel & meridian
elseif geography == 1
    PSparallel = -71; PSmeridian = 0; %Antarctic PS standard parallel & meridian
end

%create a mat-file containing geoid elevations over the ROI containing your DEMs on a 500 m-resolution grid
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
[psx_grid,psy_grid] = wgs2ps(ylon,ylat,'StandardParallel',PSparallel,'StandardMeridian',PSmeridian);

geoid.x = psx_grid; geoid.y = psy_grid; geoid.z = geoid_z;
geoid.lon = WGSlon-360; geoid.lat = WGSlat;
save([dir_output,region_abbrev,'_geoid_heights.mat'],'geoid','-v7.3');
disp('geoid heights saved');

end