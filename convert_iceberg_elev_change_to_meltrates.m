function [SL] = convert_iceberg_elev_change_to_meltrates(DEM1,DEM2,IM1,IM2,berg_numbers,geography,region_name,region_abbrev,dir_output,dir_code,dir_iceberg,dir_SMB)
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

%specify polar projection parameters
if geography == 0
    PSparallel = 70; PSmeridian = -45; %Greenland PS standard parallel & meridian
    PSprojfile = 'EPSG3413.prj'; %NSIDC Sea Ice Polar Stereographic projection proj4 text description
elseif geography == 1
    PSparallel = -71; PSmeridian = 0; %Antarctic PS standard parallel & meridian
    PSprojfile = 'EPSG3031.prj'; %Antarctic Polar Stereographic projection proj4 text description
end

%set densities (if Antarctica, iceberg densities will differ from rho_i)
rho_sw = 1026;  rho_sw_err = 2; %kg m^-3
rho_i = 900; rho_i_err = 10; %kg m^-3 

coord_files = dir([dir_iceberg,'*PScoords.txt']);
for j = 1:length(coord_files)
    coords=readmatrix(coord_files(j).name);
    PSy_early(j)=coords(1);  PSx_early(j)=coords(2); PSy_late(j)=coords(3); PSx_late(j)=coords(4);
    clear coords;
end
berg_x = nanmean([PSx_early PSx_late]); berg_y = nanmean([PSy_early PSy_late]);
%make sure DEM date string is long enough to patch now-corrected error in convert_PGC_tifs_to_matfiles.m (corrected 21/04/23)
if length(DEM1.YYYYMMDDhhmmss) < 14; DEM1.YYYYMMDDhhmmss = [DEM1.YYYYMMDDhhmmss,num2str(zeros(1,14-length(DEM1.YYYYMMDDhhmmss)))]; end 
if length(DEM2.YYYYMMDDhhmmss) < 14
    if length(DEM2.YYYYMMDDhhmmss) == 12
        DEM2.YYYYMMDDhhmmss = [DEM2.YYYYMMDDhhmmss,'00']; 
    elseif length(DEM2.YYYYMMDDhhmmss) == 13
        DEM2.YYYYMMDDhhmmss = [DEM2.YYYYMMDDhhmmss,'0']; 
    end 
end
berg_dates = [DEM1.YYYYMMDDhhmmss; DEM2.YYYYMMDDhhmmss];
to = berg_dates(1,:); tf = berg_dates(2,:);

%extract air temp (& firn density as needed) from model
density_z = [0:1:1000]; %thickness profile for density curve fitting
if geography == 1 %surface air temp, runoff, and firn density for Antarctica
    [dt,iceberg_avgtemp,surfmelt,firnair,density,f,ci] = extract_RACMO_params(dir_SMB,geography,berg_x,berg_y,berg_dates);
    density.nineseventeen = -f.b*log(-(916.9-917)/(917-f.a)); %find depth where rho=916.9 (goes to infinity at 917)
    clear FAC; FAC(1) = firnair.median; FAC(2) = firnair.median-firnair.uncert; FAC(3) = firnair.median+firnair.uncert; %estimate firn air content
    density_profile(1,:) = rho_i-(rho_i-f.a)*exp(-density_z/f.b);
    density_profile(2,:) = rho_i-(rho_i-ci(1,1))*exp(-density_z/ci(1,2)); %MINIMUM
    density_profile(3,:) = rho_i-(rho_i-ci(2,1))*exp(-density_z/ci(2,2)); %MAXIMUM
    %calculate wet density profile by flipping the shape of the exponential curve & compressing the range from (830-0) to (1026-830)
    wetdensity_profile(1,:) = 830+((830-density_profile)./830).*(rho_sw-830); wetdensity_profile(1,ceil(density.eightthir)+1:end) = density_profile(1,ceil(density.eightthir)+1:end);
    wetdensity_profile(2,:) = 830+((830-density_profile(2,:))./830).*(rho_sw-830); wetdensity_profile(2,ceil(density.eightthir)+1:end) = density_profile(2,ceil(density.eightthir)+1:end);
    wetdensity_profile(3,:) = 830+((830-density_profile(3,:))./830).*(rho_sw-830); wetdensity_profile(3,ceil(density.eightthir)+1:end) = density_profile(3,ceil(density.eightthir)+1:end);
    
    %save the FAC & density data
    if ~exist([dir_output,'firn_data/'],'dir')
        mkdir([dir_output,'firn_data/']);
    end
    save([dir_output,'firn_data/',region_name,'_density_data.mat'],'firnair','density');
    close all; drawnow;
else %only surface air temp and runoff for Greenland
    [dt,iceberg_avgtemp,surfmelt] = extract_MAR_params(dir_SMB,geography,berg_x,berg_y,berg_dates);
end

% %load the saved data if restarting
% cd(dir_iceberg);
% load([region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat']);

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
    for j = 1:length(IB)
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
    SL(i).initial.time = to; SL(i).final.time = tf; SL(i).days = dt;
    SL(i).SMB = -abs(surfmelt); SL(i).airtemp = iceberg_avgtemp;
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

%loop
for i = 1:size(SL,2)
    disp(['Extracting iceberg size info for iceberg #',num2str(i),' of ',num2str(size(SL,2))]);
    
    %zoom in
    disp(['...zooming-in & plotting the DEM ROI in the image and DEM for iceberg #',num2str(i)]);
    figure(figure2); set(gca,'ydir','normal','clim',[0 50]); hold on;
    figure(figure1); colormap gray; set(gca,'clim',[1.05*min(A.z(~isnan(A.z))) (0.95)*max(max(A.z))]);
%     imagesc(A.x,A.y,A.z); set(gca,'ydir','normal'); hold on;
    vxi = nearestneighbour(SL(i).initial.x,A.x); vyi = nearestneighbour(SL(i).initial.y,A.y);
    set(gca,'xlim',[min(A.x(vxi))-150 max(A.x(vxi))+150],'ylim',[min(A.y(vyi))-150 max(A.y(vyi))+150]);
    zoomin = questdlg('Widen the zoom window?',...
            'Zoom check','1) Yes','2) No','2) No');
        switch zoomin
            case '1) Yes'
                figure(figure1); set(gca,'xlim',[min(A.x(vxi))-250 max(A.x(vxi))+250],'ylim',[min(A.y(vyi))-250 max(A.y(vyi))+250]);
                figure(figure2); set(gca,'xlim',[min(SL(i).initial.x)-250 max(SL(i).initial.x)+250],'ylim',[min(SL(i).initial.y)-250 max(SL(i).initial.y)+250]);
            case '2) No'
                figure(figure2); set(gca,'xlim',[min(SL(i).initial.x)-150 max(SL(i).initial.x)+150],'ylim',[min(SL(i).initial.y)-150 max(SL(i).initial.y)+150]);
        end
    figure(figure1); plot(A.x(vxi),A.y(vyi),'--r','linewidth',2);
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
    
    %convert elevations over the iceberg area to a volume
    if geography == 1 %iteratively converge on best-estimate for density accounting for water saturation as needed
        %flag the iceberg as upright or overturned
        answer = questdlg('Is the iceberg upright (i.e., does it look like the glacier surface)?',...
            'Iceberg Upright or Flipped','1) Yes','2) No','1) Yes');
        switch answer
            case '1) Yes'
                SL(i).orientation = 1; SL(i).density_type = 'dry'; %upright
            case '2) No'
                SL(i).orientation = 0; SL(i).density_type = 'wet'; %overturned or a fragment of the full thickness
        end
        
        %iteratively estimate bulk density
        disp('estimating density');
        berg_densities = estimate_iceberg_density(SL(i).orientation,SL(i).initial.z_median,SL(i).initial.z_mad,density_z,density,density_profile,wetdensity_profile);
        SL(i).initial.density = berg_densities(1); %SL(i).initial.density_uncert = berg_densities_uncert(1);
        SL(i).initial_range.density = [berg_densities(2) berg_densities(3)];
        clear berg_densities;
        
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
    else
        %flag the iceberg as upright or overturned
        answer = questdlg('Is the iceberg upright (i.e., does it look like the glacier surface)?',...
            'Iceberg Upright or Flipped','1) Yes','2) No','1) Yes');
        switch answer
            case '1) Yes'
                SL(i).orientation = 1;
            case '2) No'
                SL(i).orientation = 0;
        end
        SL(i).density_type = 'N/A'; %no firn to saturate!
        
        %add density estimates to structure
        SL(i).initial.density = rho_i; %SL(i).initial.density_uncert = berg_densities_uncert(1);
        SL(i).initial_range.density = [rho_i-rho_i_err rho_i+rho_i_err];
        
        %estimate draft & thickness range based on density uncertainty
        for l = 1:2
            draft(l) = (SL(i).initial_range.density(l)/(rho_sw-SL(i).initial_range.density(l)))*SL(i).initial.z_median;
            Hberg(l) = (rho_sw/(rho_sw-SL(i).initial_range.density(l)))*SL(i).initial.z_median;
        end
    end
    SL(i).initial_range.draft = sort(draft);
    lat_area = sort(draft).*sum(dist); SL(i).initial_range.LA = lat_area;
    area = base_area*ones(size(lat_area)) + lat_area; SL(i).initial_range.TA = area;
    SL(i).initial.V = DEM_pixel_area*(rho_sw/(rho_sw-SL(i).initial.density))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(i).initial.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(i).initial_range.V(1) = DEM_pixel_area*(rho_sw/(rho_sw-min(SL(i).initial_range.density)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(i).initial.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(i).initial_range.V(2) = DEM_pixel_area*(rho_sw/(rho_sw-max(SL(i).initial_range.density)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(i).initial.z_median*sum(sum(isnan(DEM_z_masked))));
    
    %save the iceberg outline
    disp('Writing shapefile containing early iceberg polygon information');
    S.Geometry = 'Polygon';
    S.BoundingBox = [min(x) min(y); max(x) max(y)];
    S.X = double(x'); S.Y = double(y');
    S.Name = ['iceberg',num2str(berg_numbers(i).name(8:9))];
    shapefile_name = ['WV_',num2str(to),'_icebergshape',num2str(berg_numbers(i).name(8:9))];
    cd([dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);
    shapewrite(S,shapefile_name);
    copyfile([dir_code,PSprojfile],[dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_shapes/',shapefile_name,'.prj']);
    cd([dir_output,DEM1.time,'-',DEM2.time,'/']);
    clear S;

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
clear Z x1 x2 y1 y2 xlims ylims xmin xmax ymin ymax DEM_x DEM_y DEM_z ft f ci;
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
    
    %convert elevations over the iceberg area to a volume
    if geography == 1 %iteratively converge on best-estimate for density accounting for water saturation as needed
        %iteratively estimate bulk density
        disp('estimating density');
        berg_densities = estimate_iceberg_density(SL(i).orientation,SL(i).final.z_median,SL(i).final.z_mad,density_z,density,density_profile,wetdensity_profile);
        SL(i).final.density = berg_densities(1); %SL(i).initial.density_uncert = berg_densities_uncert(1);
        SL(i).final_range.density = [berg_densities(2) berg_densities(3)];
        clear berg_densities;
        
        %estimate the iceberg depth & submerged area using densities from the Ligtenberg FDM
        disp('estimating thickness & volume');
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
    else
        %add density estimates to structure
        SL(i).final.density = rho_i; %SL(i).initial.density_uncert = berg_densities_uncert(1);
        SL(i).final_range.density = [rho_i-rho_i_err rho_i+rho_i_err];
        
        %estimate draft & thickness range based on density uncertainty
        for l = 1:2
            draft(l) = (SL(i).final_range.density(l)/(rho_sw-SL(i).final_range.density(l)))*SL(i).final.z_median;
            Hberg(l) = (rho_sw/(rho_sw-SL(i).final_range.density(l)))*SL(i).final.z_median;
        end
    end
    SL(i).final_range.draft = sort(draft);
    lat_area = sort(draft).*sum(dist); SL(i).final_range.LA = lat_area;
    area = base_area*ones(size(lat_area)) + lat_area; SL(i).final_range.TA = area;
    SL(i).final.V = DEM_pixel_area*(rho_sw/(rho_sw-SL(i).final.density))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(i).final.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(i).final_range.V(1) = DEM_pixel_area*(rho_sw/(rho_sw-min(SL(i).final_range.density)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(i).final.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(i).final_range.V(2) = DEM_pixel_area*(rho_sw/(rho_sw-max(SL(i).final_range.density)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(i).final.z_median*sum(sum(isnan(DEM_z_masked))));
    
    %save the iceberg outline
    S.Geometry = 'Polygon';
    S.BoundingBox = [min(x) min(y); max(x) max(y)];
    S.X = double(x'); S.Y = double(y');
    S.Name = ['iceberg',num2str(berg_numbers(i).name(8:9))];
    shapefile_name = ['WV_',num2str(tf),'_icebergshape',num2str(berg_numbers(i).name(8:9))];
    cd([dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);
    shapewrite(S,shapefile_name);
    copyfile([dir_code,PSprojfile],[dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/',shapefile_name,'.prj']);
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
for i = 1:length(SL)
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
    if SL(i).airtemp <= 263
        rf=rf_o*exp((-60000/8.314)*((1/SL(i).airtemp)-(1/263))); %rate factor for cold ice
    else
        rf=rf_o*exp((-134000/8.314)*((1/SL(i).airtemp)-(1/263))); %rate factor for nearly-temperate ice
    end
    SL(i).ratefactor = rf;
    B = rf^(-1/3); %Pa s^1/3
    creep = ((-1/(2*sqrt(3)))*((nanmean([SL(i).initial.density SL(i).final.density])*9.81*SL(i).mean.z)/(2*sqrt(3)))^3*(1-(nanmean([SL(i).initial.density SL(i).final.density])/rho_sw))^3)/(B^3); %creep thinning rate (1/s)
%     SL(i).creep_dz = (SL(i).mean.H*creep*(dt*31536000)); %correct if dt is in decimal YEARS
    SL(i).creep_dz = (SL(i).mean.H*(creep*86400)*dt); %correct if dt is in days
    dH_submelt = dH_SMBadjust_mean + SL(i).creep_dz; %integrate creep over the ice thickness & over the time period

    %dH uncertainty sources
    %S1) Bulk density (use range in density estimates, including ice density if overturned)
    %S2) Ocean water density (use 1026 +/- 2 kg/m^3)
    %S3) Runoff uncertainty (from RACMO output)
    %R1) DEM uncertainty (random errors estimated as 2.9m)
    %R2) Translocation & Rotation (User) error (quantify by repeating procedure 10 times for each glacier)
    
    %%quantify potential bias from systematic errors
    %uncertainty from S1 & S2
    if strcmp(SL(i).density_type,'wet')
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
    for k = 1:length(IB)
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
    
end

%save the compiled data
disp('Saving melt rates');
save([dir_output,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat'],'SL','-v7.3');
disp('Melt rates saved! Advance to the next dataset');

close all; drawnow;

end

