function [SL] = update_iceberg_meltrates(DEM1,DEM2,IM1,IM2,geography,region_name,region_abbrev,iceberg_refs,dir_output,dir_iceberg,dir_code,dir_SMB)
% Function to convert iceberg elevation change to melt rates in Antarctica
% Ellyn Enderlin & Rainey Aberle, Spring 2022
%
% INPUTS:   DEM1            structure variable containing earlier DEM info
%           DEM2            structure variable containing later DEM info
%           IM1             structure variable containing earlier
%                           orthoimage info
%           IM2             structure variable containing later orthoimage
%                           info
%           iceberg_refs    number of iceberg to detect elevation change
%           dir_output      directory where all output files will be placed
%           dir_iceberg     directory where all intermediate files for
%                           individual icebergs are located
%           region_abbrev   abbrevation of region used in file names
%           region_name     name of region used in firn file name & in SL
%                           structure
%
% OUTPUTS:  Rewrites individual iceberg elevation change files (iceberg##_dz.mat) &
% updates the SL structure containing iceberg melt volume & melt rate data
close all; drawnow;

%specify polar projection parameters
if geography == 0
    PSparallel = 70; PSmeridian = -45; %Greenland PS standard parallel & meridian
    PSprojfile = 'EPSG3413.prj'; %NSIDC Sea Ice Polar Stereographic projection proj4 text description
elseif geography == 1
    PSparallel = -71; PSmeridian = 0; %Antarctic PS standard parallel & meridian
    PSprojfile = 'EPSG3031.prj'; %Antarctic Polar Stereographic projection proj4 text description
end

%specify densities
rho_sw = 1026;  rho_sw_err = 2; %kg m^-3
rho_i = 900; rho_i_err = 10; %kg m^-3

%load the saved iceberg data in the SL structure
DEM1_date = DEM1.time(1:8); DEM2_date = DEM2.time(1:8);
load([dir_output,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat']);
for i = 1:length(iceberg_refs)
    for k = 1:length(SL)
        if iceberg_refs(i) < 10
            if strmatch([num2str(0),num2str(iceberg_refs(i))],SL(k).name(end-1:end))
                berg_refs(i) = k;
            end
        else
            if strmatch(num2str(iceberg_refs(i)),SL(k).name(end-1:end))
                berg_refs(i) = k;
            end
        end
    end
end
%make sure DEM date string is long enough to patch now-corrected error in convert_PGC_tifs_to_matfiles.m (corrected 21/04/23)
if length(DEM1.YYYYMMDDhhmmss) < 14; DEM1.YYYYMMDDhhmmss = [DEM1.YYYYMMDDhhmmss,num2str(zeros(1,14-length(DEM1.YYYYMMDDhhmmss)))]; end 
if length(DEM2.YYYYMMDDhhmmss) < 14; DEM2.YYYYMMDDhhmmss = [DEM2.YYYYMMDDhhmmss,num2str(zeros(1,14-length(DEM2.YYYYMMDDhhmmss)))]; end 
berg_dates = [DEM1.YYYYMMDDhhmmss; DEM2.YYYYMMDDhhmmss];
to = berg_dates(1,:); tf = berg_dates(2,:);

%extract DEM pixel areas
DEM1_pixel_area = abs(DEM1.x(1)-DEM1.x(2)).*abs(DEM1.y(1)-DEM1.y(2)); DEM2_pixel_area = abs(DEM2.x(1)-DEM2.x(2)).*abs(DEM2.y(1)-DEM2.y(2)); %square meters

%set image bounds & crop
%EARLIER DATE
xo = []; yo = [];
for i = berg_refs
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
clear x1 x2 y1 y2 xlims ylims xmin xmax ymin ymax;
%LATER DATE
xf = []; yf = [];
for i = berg_refs
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
B.x = IM2.x(1,[xlims(1):xlims(2)]); B.y = IM2.y(1,[ylims(1):ylims(2)]);
B.z = double(IM2.z(ylims(1):ylims(2),xlims(1):xlims(2)));
if mean(B.z(~isnan(B.z)),'omitnan') < (1/3)*255 %image is dark
    z_adjust = imadjust(B.z./255,[],[],0.25);
    B.z = z_adjust; clear z_adjust;
end

%extract image pixel areas
im1_pixel_area = abs(A.x(1)-A.x(2)).*abs(A.y(1)-A.y(2)); im2_pixel_area = abs(B.x(1)-B.x(2)).*abs(B.y(1)-B.y(2)); %square meters


%update iceberg data for the icebergs for which you re-ran volume change calculations
cd(dir_iceberg);
for i = 1:length(iceberg_refs)
    if iceberg_refs(i) < 10
        berg_number = [num2str(0),num2str(iceberg_refs(i))];
    else
        berg_number = num2str(iceberg_refs(i));
    end
    berg_ref = berg_refs(i);

    %step 1: Re-extract elevation change
    disp('You should have already re-extracted elevation changes if you ran the wrapper & indicated you''re updating data');
    clear IM*; %clear full-extent images to speed up computation time
    disp('Loaded supporting data. Grabbing elevation change data.');

    %step 2: Load the re-extracted elevation change data & recalculate volume flux
    load(['iceberg',berg_number,'_dz.mat']);
    SL(berg_ref).name = [region_name,num2str(berg_number)]; 
    xo = []; yo = []; xf = []; yf = []; zo = [] zf = [];
    for j = 1:length(IB)
        xo = cat(1,xo,IB(j).vertices.xo); yo = cat(1,yo,IB(j).vertices.yo);
        xf = cat(1,xf,IB(j).vertices.xf); yf = cat(1,yf,IB(j).vertices.yf);
        non_nans = ~isnan(IB(j).zo.local_adjust.map) & ~isnan(IB(j).zf.local_adjust.rotmap);
        IB(j).zo.local_adjust.map(~non_nans) = NaN; IB(j).zf.local_adjust.rotmap(~non_nans) = NaN;
        zo(:,:,j) = IB(j).zo.local_adjust.map; zf(:,:,j) = IB(j).zf.local_adjust.rotmap;
        clear non_nans;
    end
    SL(berg_ref).initial.x = nanmean(xo,1); SL(berg_ref).initial.y = nanmean(yo,1);
    SL(berg_ref).final.x = nanmean(xf,1); SL(berg_ref).final.y = nanmean(yf,1);
    poly1 = roipoly(DEM1.x,DEM1.y,DEM1.z_elpsd_adjust,SL(berg_ref).initial.x,SL(berg_ref).initial.y);
    mask1 = logical(poly1);
    SL(berg_ref).initial.z = single(mask1.*DEM1.z_elpsd_adjust); SL(berg_ref).initial.z(SL(berg_ref).initial.z == 0) = NaN; SL(berg_ref).initial.z(SL(berg_ref).initial.z<0) = NaN;
    poly2 = roipoly(DEM2.x,DEM2.y,DEM2.z_elpsd_adjust,SL(berg_ref).final.x,SL(berg_ref).final.y);
    mask2 = logical(poly2);
    SL(berg_ref).final.z = single(mask2.*DEM2.z_elpsd_adjust); SL(berg_ref).final.z(SL(berg_ref).final.z == 0) = NaN; SL(berg_ref).final.z(SL(berg_ref).final.z<0) = NaN;
    SL(berg_ref).initial.coreg_z = IB(1).local_adjust_o; SL(berg_ref).final.coreg_z = IB(1).local_adjust_f;
    if isfield(IB,'flag')
        SL(berg_ref).flag = IB(1).flag;
    else
        SL(berg_ref).flag = NaN;
    end
    
    %plot the early DEM & image
    disp('Plotting the early image-DEM pair for the fjord');
    figure1 = figure;
    imagesc(A.x,A.y,A.z); hold on;
    colormap gray; set(gca,'clim',[1.05*min(A.z(~isnan(A.z))) (0.95)*max(max(A.z))]);
    set(gca,'ydir','normal'); hold on;
    set(gcf,'position',[50 100 600 600]);
    figure2 = figure;
    set(gcf,'position',[650 100 600 600]);
    z = DEM1.z_elpsd_adjust+SL(berg_ref).initial.coreg_z*ones(size(DEM1.z_elpsd_adjust));
    idx = nearestneighbour([A.x(1) A.x(end)],DEM1.x); idy = nearestneighbour([A.y(1) A.y(end)],DEM1.y);
    DEM_x = DEM1.x(idx(1):idx(2)); DEM_y = DEM1.y(idy(1):idy(2));
    DEM_z = z(idy(1):idy(2),idx(1):idx(2));
    imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    
    %measure the iceberg in the early image & DEM
    disp('Measure the iceberg aerial extent by drawing a polygon around the iceberg edge');
    cd([dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);
    disp(['Zooming-in & plotting the DEM ROI in the image and DEM for iceberg #',num2str(berg_ref)]);
    figure(figure1);
    %     imagesc(A.x,A.y,A.z); set(gca,'ydir','normal'); hold on;
    colormap gray; set(gca,'clim',[1.05*min(A.z(~isnan(A.z))) (0.95)*max(max(A.z))]);
    vxi = nearestneighbour(SL(berg_ref).initial.x,A.x); vyi = nearestneighbour(SL(berg_ref).initial.y,A.y);
    set(gca,'xlim',[min(A.x(vxi))-150 max(A.x(vxi))+150],'ylim',[min(A.y(vyi))-150 max(A.y(vyi))+150]);
    plot(A.x(vxi),A.y(vyi),'--r','linewidth',2);
    figure(figure2);
    %     imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    set(gca,'xlim',[min(SL(berg_ref).initial.x)-150 max(SL(berg_ref).initial.x)+150],'ylim',[min(SL(berg_ref).initial.y)-150 max(SL(berg_ref).initial.y)+150]);
    plot(SL(berg_ref).initial.x,SL(berg_ref).initial.y,'-k','linewidth',2); hold on; plot(SL(berg_ref).initial.x,SL(berg_ref).initial.y,'--w','linewidth',2); hold on;
    
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
    SL(berg_ref).initial.x_perim = x; SL(berg_ref).initial.y_perim = y;
    
    %modify the mask using the DEM if necessary
    disp('Now modify the iceberg mask needed to extract elevations from the entire iceberg area');
    clear DEM_z_masked iceberg_DEMmask x y;
    figure(figure2);
    [iceberg_DEMmask,x,y] = roipoly;
    SL(berg_ref).initial.x_perim = []; SL(berg_ref).initial.y_perim = [];
    SL(berg_ref).initial.x_perim = x; SL(berg_ref).initial.y_perim = y;
    DEM_z_masked = iceberg_DEMmask.*DEM_z;
    DEM_z_masked(iceberg_DEMmask == 0) = 0;
    clf(figure2); drawnow;
    imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    set(gca,'xlim',[min(SL(berg_ref).initial.x)-150 max(SL(berg_ref).initial.x)+150],'ylim',[min(SL(berg_ref).initial.y)-150 max(SL(berg_ref).initial.y)+150]);
    plot(SL(berg_ref).initial.x,SL(berg_ref).initial.y,'-k','linewidth',2); hold on; plot(SL(berg_ref).initial.x,SL(berg_ref).initial.y,'--w','linewidth',2); hold on;
    plot(SL(berg_ref).initial.x_perim,SL(berg_ref).initial.y_perim,'-w','linewidth',2); hold on; plot(SL(berg_ref).initial.x_perim,SL(berg_ref).initial.y_perim,'--m','linewidth',2); hold on;
    
    %use the mask to extract size & shape info
    disp('Pulling plan-view info');
    for k = 1:size(x,1)-1
        dist(k) = ((x(k)-x(k+1))^2+(y(k)-y(k+1))^2)^(1/2);
    end
    dist(end+1) = ((x(end)-x(1))^2+(y(end)-y(1))^2)^(1/2);
    shapes_o = regionprops(iceberg_IMmask,'Area','MinorAxisLength','MajorAxisLength'); %use the image-derived mask to estimate area
    lmin_o = abs(A.x(2)-A.x(1))*shapes_o.MinorAxisLength;
    lmax_o = abs(A.x(2)-A.x(1))*shapes_o.MajorAxisLength;
    SL(berg_ref).initial.width_min = lmin_o; SL(berg_ref).initial.width_max = lmax_o;
    base_area = shapes_o.Area.*im1_pixel_area; SL(berg_ref).initial.SA = base_area;
    SL(berg_ref).initial.z_map = single(DEM_z_masked); SL(berg_ref).initial.z_map(SL(berg_ref).initial.z_map == 0) = NaN; SL(berg_ref).initial.z_map(SL(berg_ref).initial.z_map<0) = NaN;
    zo = SL(berg_ref).initial.z_map(~isnan(SL(berg_ref).initial.z_map));
    SL(berg_ref).initial.z_median = nanmedian(zo); SL(berg_ref).initial.z_mad = mad(zo,1);
    perimeter = sum(dist);
    SL(berg_ref).initial.radius = perimeter/(2*pi);
    clear iceberg_IMmask;
    
    
    %extract air temp (& firn density as needed) from model
    density_z = [0:1:1000]; %thickness profile for density curve fitting
    berg_x = nanmean([SL(berg_ref).initial.x SL(berg_ref).final.x]); berg_y = nanmean([SL(berg_ref).initial.y SL(berg_ref).final.y]);
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
    SL(berg_ref).days = dt;
    SL(berg_ref).SMB = -abs(surfmelt); SL(berg_ref).airtemp = iceberg_avgtemp;
    
    %convert elevations over the iceberg area to a volume: iteratively converge on best-estimate for density accounting for water saturation as needed
    if geography == 1 
        %flag the iceberg as upright or overturned
        answer = questdlg('Is the iceberg upright (i.e., does it look like the glacier surface)?',...
            'Iceberg Upright or Flipped','1) Yes','2) No','1) Yes');
        switch answer
            case '1) Yes'
                SL(berg_ref).orientation = 1; SL(berg_ref).density_type = 'dry'; %upright
            case '2) No'
                SL(berg_ref).orientation = 0; SL(berg_ref).density_type = 'wet'; %overturned or a fragment of the full thickness
        end
        
        %iteratively estimate bulk density
        disp('estimating density');
        berg_densities = estimate_iceberg_density(SL(berg_ref).orientation,SL(berg_ref).initial.z_median,SL(berg_ref).initial.z_mad,density_z,density,density_profile,wetdensity_profile);
        SL(berg_ref).initial.density = berg_densities(1); %SL(berg_ref).initial.density_uncert = berg_densities_uncert(1);
        SL(berg_ref).initial_range.density = [berg_densities(2) berg_densities(3)];
        clear berg_densities;
        
        %estimate the iceberg depth & submerged area using densities from the Ligtenberg FDM
        disp('estimating thickness & volume');
        if SL(berg_ref).orientation == 0 %overturned/fragment
            rho_f = sort([rho_i SL(berg_ref).initial.density]); %figure out if bubble-free ice or water-saturated firn has a lower density
            for l = 1:2
                draft(l) = (rho_f(l)/(rho_sw-rho_f(l)))*SL(berg_ref).initial.z_median;
                Hberg(l) = (rho_sw/(rho_sw-rho_f(l)))*SL(berg_ref).initial.z_median;
            end
        else
            rho_f = sort(SL(berg_ref).initial_range.density); %assume density constrained by model FAC range
            for l = 1:2
                draft(l) = (rho_f(l)/(rho_sw-rho_f(l)))*SL(berg_ref).initial.z_median;
                Hberg(l) = (rho_sw/(rho_sw-rho_f(l)))*SL(berg_ref).initial.z_median;
            end
        end
    else
        %flag the iceberg as upright or overturned
        answer = questdlg('Is the iceberg upright (i.e., does it look like the glacier surface)?',...
            'Iceberg Upright or Flipped','1) Yes','2) No','1) Yes');
        switch answer
            case '1) Yes'
                SL(berg_ref).orientation = 1;
            case '2) No'
                SL(berg_ref).orientation = 0;
        end
        SL(berg_ref).density_type = 'N/A'; %no firn to saturate!
        
        %add density estimates to structure
        SL(berg_ref).initial.density = rho_i; %SL(berg_ref).initial.density_uncert = berg_densities_uncert(1);
        SL(berg_ref).initial_range.density = [rho_i-rho_i_err rho_i+rho_i_err];
        
        %estimate draft & thickness range based on density uncertainty
        for l = 1:2
            draft(l) = (SL(berg_ref).initial_range.density(l)/(rho_sw-SL(berg_ref).initial_range.density(l)))*SL(berg_ref).initial.z_median;
            Hberg(l) = (rho_sw/(rho_sw-SL(berg_ref).initial_range.density(l)))*SL(berg_ref).initial.z_median;
        end
    end
    SL(berg_ref).initial_range.draft = sort(draft);
    lat_area = sort(draft).*sum(dist); SL(berg_ref).initial_range.LA = lat_area;
    area = base_area*ones(size(lat_area)) + lat_area; SL(berg_ref).initial_range.TA = area;
    SL(berg_ref).initial.V = DEM1_pixel_area*(rho_sw/(rho_sw-SL(berg_ref).initial.density))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(berg_ref).initial.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(berg_ref).initial_range.V(1) = DEM1_pixel_area*(rho_sw/(rho_sw-min(SL(berg_ref).initial_range.density)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(berg_ref).initial.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(berg_ref).initial_range.V(2) = DEM1_pixel_area*(rho_sw/(rho_sw-max(SL(berg_ref).initial_range.density)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(berg_ref).initial.z_median*sum(sum(isnan(DEM_z_masked))));
    
    
    %save the iceberg outline
    disp('Writing shapefile containing early iceberg polygon information');
    S.Geometry = 'Polygon';
    S.BoundingBox = [min(x) min(y); max(x) max(y)];
    S.X = double(x'); S.Y = double(y');
    S.Name = ['iceberg',num2str(berg_number)];
    shapefile_name = ['WV_',num2str(to),'_icebergshape',num2str(berg_number)];
    cd([dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);
    shapewrite(S,shapefile_name);
    copyfile([dir_code,PSprojfile],[dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_shapes/',shapefile_name,'.prj']);
    cd([dir_output,DEM1.time,'-',DEM2.time,'/']);
    clear S;
    
    %save the data
    disp('Saving data');
    save([dir_output,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat'],'SL','-v7.3');
    cd([dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);
    clear DEM_z_masked z idx idy x y iceberg_IMmask iceberg_DEMmask rho_f draft Hberg lat_area base_area dist area_mask masked_pixels image_xgrid image_ygrid shapes lmin* lmax* perimeter;
    clear x1 x2 y1 y2 xlims ylims xmin xmax ymin ymax DEM_x DEM_y DEM_z;
    disp('Advancing to the later date');
    close all; drawnow;
    
    %plot the later DEM & image
    disp('Plotting the later image-DEM pair for the fjord');
    figure1 = figure;
    imagesc(B.x,B.y,B.z); hold on;
    colormap gray; set(gca,'clim',[1.05*min(B.z(~isnan(B.z))) (0.95)*max(max(B.z))]);
    set(gca,'ydir','normal'); hold on;
    set(gcf,'position',[50 100 600 600]);
    figure2 = figure;
    set(gcf,'position',[650 100 600 600]);
    z = DEM2.z_elpsd_adjust+SL(berg_ref).final.coreg_z*ones(size(DEM2.z_elpsd_adjust));
    idx = nearestneighbour([B.x(1) B.x(end)],DEM2.x); idy = nearestneighbour([B.y(1) B.y(end)],DEM2.y);
    DEM_x = DEM2.x(idx(1):idx(2)); DEM_y = DEM2.y(idy(1):idy(2));
    DEM_z = z(idy(1):idy(2),idx(1):idx(2));
    imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    
    %measure the iceberg in the later image & DEM
    cd([dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);
    disp(['Zooming-in & plotting the DEM ROI in the image and DEM for iceberg #',num2str(berg_ref)]);
    figure(figure1);
    colormap gray; set(gca,'clim',[1.05*min(B.z(~isnan(B.z))) (0.95)*max(max(B.z))]);
    vxf = nearestneighbour(SL(berg_ref).final.x,B.x); vyf = nearestneighbour(SL(berg_ref).final.y,B.y);
    set(gca,'xlim',[min(SL(berg_ref).final.x)-150 max(SL(berg_ref).final.x)+150],'ylim',[min(SL(berg_ref).final.y)-150 max(SL(berg_ref).final.y)+150]);
    plot(B.x(vxf),B.y(vyf),'--r','linewidth',2);
    figure(figure2);
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    set(gca,'xlim',[min(SL(berg_ref).final.x)-150 max(SL(berg_ref).final.x)+150],'ylim',[min(SL(berg_ref).final.y)-150 max(SL(berg_ref).final.y)+150]);
    plot(SL(berg_ref).final.x,SL(berg_ref).final.y,'-k','linewidth',2); hold on; plot(SL(berg_ref).final.x,SL(berg_ref).final.y,'--w','linewidth',2); hold on;
    
    %define the iceberg area & save its vertices to a shapefile
    disp('Draw a polygon following the iceberg edges in the WV image (used to estimate surface area).');
    figure(figure1);
    [iceberg_IMmask,x,y] = roipoly;
    
    %apply the mask from the image to the DEM
    [image_xgrid,image_ygrid] = meshgrid(B.x,B.y);
    [DEM_xgrid,DEM_ygrid] = meshgrid(DEM_x,DEM_y);
    iceberg_DEMmask = interp2(image_xgrid,image_ygrid,double(iceberg_IMmask),DEM_xgrid,DEM_ygrid);
    iceberg_DEMmask = round(iceberg_DEMmask);
    DEM_z_masked = iceberg_DEMmask.*DEM_z;
    DEM_z_masked(iceberg_DEMmask == 0) = 0;
    disp('NOTE: DEM coordinates may be shifted relative to the image: compare the DEM ROI (black & white dash) & the iceberg polygon (pink & white dash)on the DEM');
    figure(figure2);
    plot(x,y,'-w','linewidth',2); hold on; plot(x,y,'--m','linewidth',2); hold on;
    SL(berg_ref).final.x_perim = x; SL(berg_ref).final.y_perim = y;
    
    %modify the mask using the DEM if necessary
    disp('Now modify the iceberg mask need to extract elevations from the entire iceberg area');
    clear DEM_z_masked iceberg_DEMmask x y;
    figure(figure2);
    [iceberg_DEMmask,x,y] = roipoly;
    SL(berg_ref).final.x_perim = []; SL(berg_ref).final.y_perim = [];
    SL(berg_ref).final.x_perim = x; SL(berg_ref).final.y_perim = y;
    DEM_z_masked = iceberg_DEMmask.*DEM_z;
    DEM_z_masked(iceberg_DEMmask == 0) = 0;
    clf(figure2); drawnow;
    imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    set(gca,'xlim',[min(SL(berg_ref).final.x)-150 max(SL(berg_ref).final.x)+150],'ylim',[min(SL(berg_ref).final.y)-150 max(SL(berg_ref).final.y)+150]);
    plot(SL(berg_ref).final.x,SL(berg_ref).final.y,'-k','linewidth',2); hold on; plot(SL(berg_ref).final.x,SL(berg_ref).final.y,'--w','linewidth',2); hold on;
    plot(SL(berg_ref).final.x_perim,SL(berg_ref).final.y_perim,'-w','linewidth',2); hold on; plot(SL(berg_ref).final.x_perim,SL(berg_ref).final.y_perim,'--m','linewidth',2); hold on;
    
    %use the mask to extract size & shape info
    disp('Pulling plan-view info');
    for k = 1:size(x,1)-1
        dist(k) = ((x(k)-x(k+1))^2+(y(k)-y(k+1))^2)^(1/2);
    end
    dist(end+1) = ((x(end)-x(1))^2+(y(end)-y(1))^2)^(1/2);
    shapes_f = regionprops(iceberg_IMmask,'Area','MinorAxisLength','MajorAxisLength');
    lmin_f = abs(B.x(2)-B.x(1))*shapes_f.MinorAxisLength;
    lmax_f = abs(B.x(2)-B.x(1))*shapes_f.MajorAxisLength;
    SL(berg_ref).final.width_min = lmin_f; SL(berg_ref).final.width_max = lmax_f;
    base_area = shapes_f.Area.*im2_pixel_area; SL(berg_ref).final.SA = base_area;
    SL(berg_ref).final.z_map = single(DEM_z_masked); SL(berg_ref).final.z_map(SL(berg_ref).final.z_map == 0) = NaN; SL(berg_ref).final.z_map(SL(berg_ref).final.z_map<0) = NaN;
    zf = SL(berg_ref).final.z_map(~isnan(SL(berg_ref).final.z_map));
    SL(berg_ref).final.z_median = nanmedian(zf); SL(berg_ref).final.z_mad = mad(zf,1);
    perimeter = sum(dist);
    SL(berg_ref).final.radius = perimeter/(2*pi);
    clear iceberg_IMmask;
    
    %convert elevations over the iceberg area to a volume
    if geography == 1 %iteratively converge on best-estimate for density accounting for water saturation as needed
        %iteratively estimate bulk density
        disp('estimating density');
        berg_densities = estimate_iceberg_density(SL(berg_ref).orientation,SL(berg_ref).final.z_median,SL(berg_ref).final.z_mad,density_z,density,density_profile,wetdensity_profile);
        SL(berg_ref).final.density = berg_densities(1); %SL(berg_ref).initial.density_uncert = berg_densities_uncert(1);
        SL(berg_ref).final_range.density = [berg_densities(2) berg_densities(3)];
        clear berg_densities;
        
        %estimate the iceberg depth & submerged area using densities from the Ligtenberg FDM
        disp('estimating thickness & volume');
        if SL(berg_ref).orientation == 0 %overturned/fragment
            rho_f = sort([rho_i SL(berg_ref).final.density]); %figure out if bubble-free ice or water-saturated firn has a lower density
            for l = 1:2
                draft(l) = (rho_f(l)/(rho_sw-rho_f(l)))*SL(berg_ref).final.z_median;
                Hberg(l) = (rho_sw/(rho_sw-rho_f(l)))*SL(berg_ref).final.z_median;
            end
        else
            rho_f = sort(SL(berg_ref).final_range.density); %assume density constrained by model FAC range
            for l = 1:2
                draft(l) = (rho_f(l)/(rho_sw-rho_f(l)))*SL(berg_ref).final.z_median;
                Hberg(l) = (rho_sw/(rho_sw-rho_f(l)))*SL(berg_ref).final.z_median;
            end
        end
    else
        %add density estimates to structure
        SL(berg_ref).final.density = rho_i; %SL(berg_ref).initial.density_uncert = berg_densities_uncert(1);
        SL(berg_ref).final_range.density = [rho_i-rho_i_err rho_i+rho_i_err];
        
        %estimate draft & thickness range based on density uncertainty
        for l = 1:2
            draft(l) = (SL(berg_ref).final_range.density(l)/(rho_sw-SL(berg_ref).final_range.density(l)))*SL(berg_ref).final.z_median;
            Hberg(l) = (rho_sw/(rho_sw-SL(berg_ref).final_range.density(l)))*SL(berg_ref).final.z_median;
        end
    end
    SL(berg_ref).final_range.draft = sort(draft);
    lat_area = sort(draft).*sum(dist); SL(berg_ref).final_range.LA = lat_area;
    area = base_area*ones(size(lat_area)) + lat_area; SL(berg_ref).final_range.TA = area;
    SL(berg_ref).final.V = DEM2_pixel_area*(rho_sw/(rho_sw-SL(berg_ref).final.density))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(berg_ref).final.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(berg_ref).final_range.V(1) = DEM2_pixel_area*(rho_sw/(rho_sw-min(SL(berg_ref).final_range.density)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(berg_ref).final.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(berg_ref).final_range.V(2) = DEM2_pixel_area*(rho_sw/(rho_sw-max(SL(berg_ref).final_range.density)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(berg_ref).final.z_median*sum(sum(isnan(DEM_z_masked))));
    
    %save the iceberg outline
    S.Geometry = 'Polygon';
    S.BoundingBox = [min(x) min(y); max(x) max(y)];
    S.X = double(x'); S.Y = double(y');
    S.Name = ['iceberg',num2str(berg_number)];
    shapefile_name = ['WV_',num2str(tf),'_icebergshape',num2str(berg_number)];
    cd([dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);
    shapewrite(S,shapefile_name);
    copyfile([dir_code,PSprojfile],[dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/',shapefile_name,'.prj']);
    cd([dir_output,DEM1.time,'-',DEM2.time,'/']);
    clear S;
    
    %calculate mean surface area for the exposed face (SA) &
    %submerged (TA) portions of the iceberg & associated temporal uncertainty
    SL(berg_ref).mean.z = (SL(berg_ref).initial.z_median+SL(berg_ref).final.z_median)/2;
    SL(berg_ref).mean.SA = (SL(berg_ref).initial.SA+SL(berg_ref).final.SA)/2;
    SL(berg_ref).change.SA = abs(SL(berg_ref).initial.SA-SL(berg_ref).final.SA)/2;
    SL(berg_ref).mean.TA = (nanmean(SL(berg_ref).initial_range.TA)+nanmean(SL(berg_ref).final_range.TA))/2;
    SL(berg_ref).uncert.TA = max([(SL(berg_ref).initial_range.TA(2)-SL(berg_ref).initial_range.TA(1))/2 (SL(berg_ref).final_range.TA(2)-SL(berg_ref).final_range.TA(1))/2]);
    SL(berg_ref).change.TA = abs(nanmean(SL(berg_ref).initial_range.TA)-nanmean(SL(berg_ref).final_range.TA))/2;
    SL(berg_ref).mean.V = (nanmean(SL(berg_ref).initial_range.V)+nanmean(SL(berg_ref).final_range.V))/2;
    SL(berg_ref).uncert.V = max([(SL(berg_ref).initial_range.V(2)-SL(berg_ref).initial_range.V(1))/2 (SL(berg_ref).final_range.V(2)-SL(berg_ref).final_range.V(1))/2]);
    SL(berg_ref).change.V = abs(nanmean(SL(berg_ref).initial_range.V)-nanmean(SL(berg_ref).final_range.V))/2;
    
    %save the updated size data
    disp('Saving data');
    save([dir_output,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat'],'SL','-v7.3');
    clear DEM_z_masked z idx idy x y iceberg_IMmask iceberg_DEMmask rho_f draft Hberg lat_area base_area dist area_mask masked_pixels image_xgrid image_ygrid shapes lmin* lmax* perimeter;
    clear x1 x2 y1 y2 xlims ylims xmin xmax ymin ymax DEM_x DEM_y DEM_z;
    close all; drawnow;
%     
    %estimate the density
    if SL(berg_ref).orientation == 0
        rho_f = nanmean([SL(berg_ref).initial.density SL(berg_ref).final.density]);
        rho_f_range = [rho_i SL(berg_ref).initial_range.density SL(berg_ref).final_range.density];
        rho_f_err = [min(rho_f_range-rho_f) max(rho_f_range-rho_f)];
    else
        rho_f = nanmean([SL(berg_ref).initial.density SL(berg_ref).final.density]);
        rho_f_range = [SL(berg_ref).initial_range.density SL(berg_ref).final_range.density];
        rho_f_err = [min(rho_f_range-rho_f) max(rho_f_range-rho_f)];
    end
    
    %calculate the ice thickness (change to H = 3*(V/SA) if assuming the submerged shape is a cone!)
    dt=sum(SL(berg_ref).days);
    SL(berg_ref).initial_range.H = (repmat(rho_sw,size(SL(berg_ref).initial_range.density))./(repmat(rho_sw,size(SL(berg_ref).initial_range.density))-SL(berg_ref).initial_range.density)).*SL(berg_ref).initial.z_median;
    SL(berg_ref).final_range.H = (repmat(rho_sw,size(SL(berg_ref).final_range.density))./(repmat(rho_sw,size(SL(berg_ref).final_range.density))-SL(berg_ref).final_range.density)).*SL(berg_ref).final.z_median;
    SL(berg_ref).mean.H = nanmean([SL(berg_ref).initial_range.H SL(berg_ref).final_range.H]);
    SL(berg_ref).change.H = nanmean(SL(berg_ref).initial_range.H) - nanmean(SL(berg_ref).final_range.H);
    SL(berg_ref).mean.draft = nanmean([SL(berg_ref).initial_range.draft SL(berg_ref).final_range.draft]);
    SL(berg_ref).change.draft = nanmean(SL(berg_ref).initial_range.draft) - nanmean(SL(berg_ref).final_range.draft);
    
    %convert surface elevation changes to thickness changes
    vals = [];
    for j = 1:length(IB)
%         if IB(j).local_adjust_f ~= 0 %use local sea level-adjusted elevations
            non_NaNs = length(IB(j).dz.local_adjust.map(~isnan(IB(j).dz.local_adjust.map)));
            vals(size(vals,2)+1:size(vals,2)+non_NaNs) = IB(j).dz.local_adjust.map(~isnan(IB(j).dz.local_adjust.map));
%         else %use bedrock & tide-adjusted elevations
%             non_NaNs = length(IB(j).dz.br_tide_adjust.map(~isnan(IB(j).dz.br_tide_adjust.map)));
%             vals(size(vals,2)+1:size(vals,2)+non_NaNs) = IB(j).dz.br_tide_adjust.map(~isnan(IB(j).dz.br_tide_adjust.map));
%          end
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
    SL(berg_ref).mean.dz = nanmean(vals); SL(berg_ref).uncert.dz = nanstd(vals);
    dZ_mean = (rho_sw/(rho_sw-nanmean([SL(berg_ref).initial.density SL(berg_ref).final.density])))*SL(berg_ref).mean.dz;
    dZ_stdev = (rho_sw/(rho_sw-nanmean([SL(berg_ref).initial.density SL(berg_ref).final.density])))*SL(berg_ref).uncert.dz;
    
    %correct for SMB loss
    dH_SMBadjust_mean = dZ_mean + SL(berg_ref).SMB; %SL(berg_ref).SMB=-surfmelt
    
    %correct for creep thinning
    rf_o = 3.5e-25; %rate factor at -10C (Pa^-3 s^-1) from Cuffey & Paterson p. 74
    if iceberg_avgtemp <= 263
        rf=rf_o*exp((-60000/8.314)*((1/SL(berg_ref).airtemp)-(1/263))); %rate factor for cold ice
    else
        rf=rf_o*exp((-134000/8.314)*((1/SL(berg_ref).airtemp)-(1/263))); %rate factor for nearly-temperate ice
    end
    SL(berg_ref).ratefactor = rf;
    B = rf^(-1/3); %Pa s^1/3
    creep = ((-1/(2*sqrt(3)))*((nanmean([SL(berg_ref).initial.density SL(berg_ref).final.density])*9.81*SL(berg_ref).mean.z)/(2*sqrt(3)))^3*(1-(nanmean([SL(berg_ref).initial.density SL(berg_ref).final.density])/rho_sw))^3)/(B^3); %creep thinning rate (1/s)
    SL(berg_ref).creep_dz = (SL(berg_ref).mean.H*(creep*86400)*dt);
    dH_submelt = dH_SMBadjust_mean + SL(berg_ref).creep_dz; %integrate creep over the ice thickness & over the time period
    
    %dH uncertainty sources
    %S1) Bulk density (use range in density estimates, including ice density if overturned)
    %S2) Ocean water density (use 1026 +/- 2 kg/m^3)
    %S3) Runoff uncertainty (from RACMO output)
    %R1) DEM uncertainty (random errors estimated as 2.9m)
    %R2) Translocation & Rotation (User) error (quantify by repeating procedure 10 times for each glacier)
    
    %%quantify potential bias from systematic errors
    %uncertainty from S1 & S2
    if strcmp(SL(berg_ref).density_type,'wet')
        rho_f = nanmean([SL(berg_ref).initial.density SL(berg_ref).final.density]);
        rho_f_range = [rho_i SL(berg_ref).initial_range.density SL(berg_ref).final_range.density];
        rho_f_err = [min(rho_f_range-rho_f) max(rho_f_range-rho_f)];
    else
        rho_f = nanmean([SL(berg_ref).initial.density SL(berg_ref).final.density]);
        rho_f_range = [SL(berg_ref).initial_range.density SL(berg_ref).final_range.density];
        rho_f_err = [min(rho_f_range-rho_f) max(rho_f_range-rho_f)];
    end
    rho_added_err = sqrt(rho_f_err.^2 + repmat(rho_sw_err,size(rho_f_err)).^2);
    rho_conversion = rho_sw/(rho_sw-rho_f);
    rho_conversion_err = abs(rho_conversion)*sqrt(repmat((rho_sw_err/rho_sw),size(rho_f_err)).^2 + (rho_added_err./repmat(rho_sw-rho_f,size(rho_f_err))).^2);
    
    %uncertainty from S3
    dz_SMB_err = 0.3*abs(SL(berg_ref).SMB);
    
    %uncertainty from R1
    pixelno_o = SL(berg_ref).initial.SA./DEM1_pixel_area; pixelno_f = SL(berg_ref).final.SA./DEM2_pixel_area;
    dz_randerr = sqrt((2.9^2+2.9^2)./nanmean([pixelno_o pixelno_f])); %(2.9^2+2.9^2) is to account for DEM differencing over time
    
    %uncertainty from R2
    for k = 1:length(IB)
        dzs(k) = IB(k).dz.local_adjust.mean;
    end
    dz_user_err = nanstd(dzs);
    
    %total, systematic, and random errors
    SL(berg_ref).uncert.h_SMB = dz_SMB_err;
    SL(berg_ref).uncert.h_rand = dz_randerr;
    SL(berg_ref).uncert.h_user = dz_user_err;
    SL(berg_ref).uncert.hH_conversion = rho_conversion_err;
    dH_total_err = abs(dH_submelt).*sqrt(repmat(((dz_SMB_err^2 + dz_randerr^2 + dz_user_err^2)/SL(berg_ref).mean.dz^2),size(rho_conversion_err)) + (rho_conversion_err./repmat(rho_conversion,size(rho_conversion_err))).^2);
    
    %convert all 'dHs', which are change in thickness as estimated using the elevation change,
    %(& errors) to changes in volume so that the iceberg melt rate
    %(dH/dt) can be calculated by dividing the volume change by the total
    %submerged iceberg area
    dV_mean = SL(berg_ref).mean.SA*dH_submelt;
    dV_stdev = abs(dV_mean).*sqrt((dZ_stdev./dH_submelt)^2);
    dV_total_err = abs(dV_mean).*sqrt((dH_total_err./repmat(dH_submelt,size(dH_total_err))).^2 + repmat((SL(berg_ref).change.SA/2)./SL(berg_ref).mean.SA,size(dH_total_err)).^2);
    
    %convert to a melt rate (volume change with time & thickness change
    %over the entire submerged face)
    dVdt_mean = dV_mean/dt;  dVdt_stdev = dV_stdev/dt;
    dVdt_total_err = dV_total_err./dt;
    dHdt_mean = dVdt_mean/SL(berg_ref).mean.TA;
    dHdt_stdev = abs(dHdt_mean).*sqrt((dVdt_stdev./dVdt_mean)^2+(SL(berg_ref).uncert.TA./SL(berg_ref).mean.TA)^2);
    dHdt_total_err = abs(dHdt_mean).*sqrt((dVdt_total_err./repmat(dVdt_mean,size(dVdt_total_err))).^2 + repmat(SL(berg_ref).uncert.TA./SL(berg_ref).mean.TA,size(dVdt_total_err)).^2);
    
    %save volume fluxes and melt rates to the structure
    SL(berg_ref).mean.dVdt = dVdt_mean;
    SL(berg_ref).range.dVdt(1) = dVdt_mean-3*dVdt_stdev; SL(berg_ref).range.dVdt(2) = dVdt_mean+3*dVdt_stdev;
    SL(berg_ref).uncert.dVdt = dVdt_total_err;
    SL(berg_ref).mean.dHdt = dHdt_mean;
    SL(berg_ref).range.dHdt(1) = dHdt_mean-3*dHdt_stdev; SL(berg_ref).range.dHdt(2) = dHdt_mean+3*dHdt_stdev;
    SL(berg_ref).uncert.dHdt = dHdt_total_err;
    
    %if overlapping bedrock regions were present, add the uniform offset that
    %would be applied from bedrock & tidal differences to the structure
%     if ~isnan(IB(1).dz.br_tide_adjust.mean)
%         SL(berg_ref).elpsd_adjust_o = IB(1).elpsd_adjust_o; SL(berg_ref).elpsd_bt_adjust_f = IB(1).elpsd_adjust_f+IB(1).dbedrock_f+IB(1).dtide_f;
%   end
    
end

%save the compiled data
disp('Saving melt rates');
save([dir_output,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat'],'SL','-v7.3');
disp('Melt rates resaved!');
close all; drawnow;
clear A B Y Z;

%plot the data to determine if any icebergs have bad melt estimates
plot_flag = 1; %plot data
table_flag = 0; %suppress export to table
plot_export_iceberg_melt_data(SL,dir_output,region_abbrev,DEM1,DEM2,plot_flag,table_flag);


end