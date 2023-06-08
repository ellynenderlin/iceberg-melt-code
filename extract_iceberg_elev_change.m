function [IB,dz] = extract_iceberg_elev_change(DEM1,DEM2,IM1,IM2,iceberg_no,dir_output,dir_code,geography,region_abbrev)
% Function to estimate iceberg freshwater fluxes and melt rates
% Ellyn Enderlin & Rainey Aberle, Fall 2021
%
% INPUTS:   DEM1            structure variable containing earlier DEM info
%           DEM2            structure variable containing later DEM info
%           IM1             structure variable containing earlier
%                           orthoimage info
%           IM2             structure variable containing later orthoimage
%                           info
%           iceberg_no      number of iceberg to detect elevation change
%           dir_output      directory where output files will be saved
%           region_abbrev   abbrevation of region used in file names
%
% OUTPUTS:  IB              structure variable containing info for each
%                               iceberg
%           dz              structure variable containing change in 
%                               elevation info
% 
% Calls the following external functions:
%   - measure_iceberg_motion.m
%   - unrotate_untranslate_iceberg.m
%   - sea_level_adjust.m
%   - nearestneighbour.m
elev_cmap = cmocean('thermal',10001); elev_cmap(1,:) = [1 1 1]; 

close all; drawnow;
%specify polar projection parameters
if geography == 0
    PSparallel = 70; PSmeridian = -45; %Greenland PS standard parallel & meridian
    PSprojfile = 'EPSG3413.prj'; %NSIDC Sea Ice Polar Stereographic projection proj4 text description
elseif geography == 1
    PSparallel = -71; PSmeridian = 0; %Antarctic PS standard parallel & meridian
    PSprojfile = 'EPSG3031.prj'; %Antarctic Polar Stereographic projection proj4 text description
end
disp(['Extracting changes in elevation for iceberg ',iceberg_no]);

%figure out the size of the monitor(s) used for plotting in order to
%determine the appropriate figure sizes and positions
monitor_sizes = get(0,'MonitorPositions');
ext_monitor = find(monitor_sizes(:,1) == max(monitor_sizes(:,1)));
if size(monitor_sizes,1) > 1
    DEMo_pos = [monitor_sizes(ext_monitor,1)+50 (5/6)*(monitor_sizes(ext_monitor,4)-50)/2+50 (5/6)*(monitor_sizes(ext_monitor,3)-50)/2 (5/6)*(monitor_sizes(ext_monitor,4)-50)/2]; 
    DEMf_pos = [monitor_sizes(ext_monitor,1)+50 0 (5/6)*(monitor_sizes(ext_monitor,3)-50)/2 (5/6)*(monitor_sizes(ext_monitor,4)-50)/2];
    imo_pos = [(monitor_sizes(ext_monitor,3)-50)/2+monitor_sizes(ext_monitor,1)+50 (5/6)*(monitor_sizes(ext_monitor,4)-50)/2+50 (5/6)*(monitor_sizes(ext_monitor,3)-50)/2 (5/6)*(monitor_sizes(ext_monitor,4)-50)/2]; 
    imf_pos = [(monitor_sizes(ext_monitor,3)-50)/2+monitor_sizes(ext_monitor,1)+50 0 (5/6)*(monitor_sizes(ext_monitor,3)-50)/2 (5/6)*(monitor_sizes(ext_monitor,4)-50)/2];
else
    DEMo_pos = [50 (5/6)*(monitor_sizes(4)-50)/2+50 (5/6)*(monitor_sizes(3)-50)/2 (5/6)*(monitor_sizes(4)-50)/2]; 
    DEMf_pos = [50 50 (5/6)*(monitor_sizes(3)-50)/2 (5/6)*(monitor_sizes(4)-50)/2];
    imo_pos = [(monitor_sizes(3)-50)/2+50 (5/6)*(monitor_sizes(4)-50)/2+50 (5/6)*(monitor_sizes(3)-50)/2 (5/6)*(monitor_sizes(4)-50)/2]; 
    imf_pos = [(monitor_sizes(3)-50)/2+50 50 (5/6)*(monitor_sizes(3)-50)/2 (5/6)*(monitor_sizes(4)-50)/2];
end

%calculate the time change between images %used to reference date variable in these structures (changed 26/05/2022)
date_o = [DEM1.time(1:4),'/',DEM1.time(5:6),'/',DEM1.time(7:8)]; 
date_f = [DEM2.time(1:4),'/',DEM2.time(5:6),'/',DEM2.time(7:8)]; 

%load the sea level offset file
fjord = load([dir_output,'fjord_offset_',DEM1.time,'-',DEM2.time,'.mat']).fjord;

%crop the data domains
DEM1_dx = abs(DEM1.x(1)-DEM1.x(2)); DEM2_dx = abs(DEM2.x(1)-DEM2.x(2));
I1_dx = abs(IM1.x(1)-IM1.x(2)); I2_dx = abs(IM2.x(1)-IM2.x(2));
DEMhalo1 = round(3000/DEM1_dx); DEMhalo2 = round(3000/DEM2_dx);
IMhalo1 = round(3000/I1_dx); IMhalo2 = round(3000/I2_dx);
icebergs_dz = dir([dir_output,'/',DEM1.time,'-',DEM2.time,'/','iceberg*dz.mat']);
PScoord_files = dir([dir_output,'/',DEM1.time,'-',DEM2.time,'/','iceberg*PScoords.txt']);
if ~isempty(PScoord_files)
    %read the coordinates for the iceberg
    iceberg_coords = [dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg',iceberg_no,'_PScoords.txt'];
    coords = cell2mat(textscan(fopen(iceberg_coords),'%f64 %f64 %f64 %f64','Delimiter',',','headerlines',1));
    early_y = coords(1); early_x = coords(2);
    late_y = coords(3); late_x = coords(4);
    
    %find the location of the iceberg center in the DEMs
    xce = nearestneighbour(early_x,DEM1.x); yce = nearestneighbour(early_y,DEM1.y);
    xc = nearestneighbour(late_x,DEM2.x); yc = nearestneighbour(late_y,DEM2.y);
    xmine = xce-DEMhalo1; xmaxe = xce+DEMhalo1; ymine = yce-DEMhalo1; ymaxe = yce+DEMhalo1;
    xmine(xmine<=1) = 1; xmaxe(xmaxe >= length(DEM1.x)) = length(DEM1.x);
    ymine(ymine<=1) = 1; ymaxe(ymaxe >= length(DEM1.y)) = length(DEM1.y);
    xmin = xc-DEMhalo2; xmax = xc+DEMhalo2; ymin = yc-DEMhalo2; ymax = yc+DEMhalo2;
    xmin(xmin<=1) = 1; xmax(xmax >= length(DEM2.x)) = length(DEM2.x);
    ymin(ymin<=1) = 1; ymax(ymax >= length(DEM2.y)) = length(DEM2.y);
    
    %find the location of the iceberg center in the images
    xcie = nearestneighbour(early_x,IM1.x); ycie = nearestneighbour(early_y,IM1.y);
    xminie = xcie-IMhalo1; xmaxie = xcie+IMhalo1; yminie = ycie-IMhalo1; ymaxie = ycie+IMhalo1;
    xminie(xminie<=1) = 1; xmaxie(xmaxie >= length(IM1.x)) = length(IM1.x);
    yminie(yminie<=1) = 1; ymaxie(ymaxie >= length(IM1.y)) = length(IM1.y);
    xci = nearestneighbour(late_x,IM2.x); yci = nearestneighbour(late_y,IM2.y);
    xmini = xci-IMhalo2; xmaxi = xci+IMhalo2; ymini = yci-IMhalo2; ymaxi = yci+IMhalo2;
    xmini(xmini<=1) = 1; xmaxi(xmaxi >= length(IM2.x)) = length(IM2.x);
    ymini(ymini<=1) = 1; ymaxi(ymaxi >= length(IM2.y)) = length(IM2.y);

    %plot the iceberg images & assess whether the center coordinates need to be adjusted
    figureA = figure; set(figureA,'position',imo_pos); %set(figureA,'position',[50 250 800 700]);
    imagesc(IM1.x(1,xminie:xmaxie),IM1.y(1,yminie:ymaxie),IM1.z(yminie:ymaxie,xminie:xmaxie)); axis xy equal; colormap gray; hold on;
    % plot iceberg coordinates on early image
    hold on; plot(early_x,early_y,'*w','markersize',7,'linewidth',1);
    set(gca,'fontsize',14);
    title(['Early date: ',datestr(date_o,'yyyy/mm/dd')],'fontsize',16);
    figureB = figure; set(figureB,'position',imf_pos); %set(figureB,'position',[850 250 800 700]);
    imagesc(IM2.x(1,xmini:xmaxi),IM2.y(1,ymini:ymaxi),IM2.z(ymini:ymaxi,xmini:xmaxi)); axis xy equal; colormap gray; hold on;
    % plot iceberg coordinates on later image
    hold on; plot(late_x,late_y,'*w','markersize',7,'linewidth',1);
    set(gca,'fontsize',14);
    title(['Late date: ',datestr(date_f,'yyyy/mm/dd')],'fontsize',16);
    
    %adjust coordinates if necessary
    disp('The iceberg of interest should be in the middle of the figures...');
    prompt = 'Do the iceberg coordinates need to be modified (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        disp('Type new values for early_y, early_x, late_y, & late_x in the command window separated by semicolons & followed by "dbcont" (w/o quotes) then hit enter to resume');
        keyboard
        coords = [early_y early_x late_y late_x];
        dlmwrite(['iceberg',iceberg_no,'_PScoords.txt'],coords,'\t');
        
        %replot images
        close all;
        xce = nearestneighbour(early_x,DEM1.x); yce = nearestneighbour(early_y,DEM1.y);
        xc = nearestneighbour(late_x,DEM2.x); yc = nearestneighbour(late_y,DEM2.y);
        xmine = xce-DEMhalo1; xmaxe = xce+DEMhalo1; ymine = yce-DEMhalo1; ymaxe = yce+DEMhalo1;
        xmine(xmine<=1) = 1; xmaxe(xmaxe >= length(DEM1.x)) = length(DEM1.x);
        ymine(ymine<=1) = 1; ymaxe(ymaxe >= length(DEM1.y)) = length(DEM1.y);
        xmin = xc-DEMhalo2; xmax = xc+DEMhalo2; ymin = yc-DEMhalo2; ymax = yc+DEMhalo2;
        xmin(xmin<=1) = 1; xmax(xmax >= length(DEM2.x)) = length(DEM2.x);
        ymin(ymin<=1) = 1; ymax(ymax >= length(DEM2.y)) = length(DEM2.y);
        
        %find the location of the iceberg center in the images
        xcie = nearestneighbour(early_x,IM1.x); ycie = nearestneighbour(early_y,IM1.y);
        xminie = xcie-IMhalo1; xmaxie = xcie+IMhalo1; yminie = ycie-IMhalo1; ymaxie = ycie+IMhalo1;
        xminie(xminie<=1) = 1; xmaxie(xmaxie >= length(IM1.x)) = length(IM1.x);
        yminie(yminie<=1) = 1; ymaxie(ymaxie >= length(IM1.y)) = length(IM1.y);
        xci = nearestneighbour(late_x,IM2.x); yci = nearestneighbour(late_y,IM2.y);
        xmini = xci-IMhalo2; xmaxi = xci+IMhalo2; ymini = yci-IMhalo2; ymaxi = yci+IMhalo2;
        xmini(xmini<=1) = 1; xmaxi(xmaxi >= length(IM2.x)) = length(IM2.x);
        ymini(ymini<=1) = 1; ymaxi(ymaxi >= length(IM2.y)) = length(IM2.y);
        
        %plot the iceberg images & assess whether the center coordinates need to be adjusted
        figureA = figure; set(figureA,'position',imo_pos); %set(figureA,'position',[50 250 800 700]);
        imagesc(IM1.x(1,xminie:xmaxie),IM1.y(1,yminie:ymaxie),IM1.z_masked(yminie:ymaxie,xminie:xmaxie)); axis xy equal; colormap gray; hold on;
        set(gca,'fontsize',14);
        title(['Early date: ',num2str(datestr(date_o,'yyyy/mm/dd'))],'fontsize',16);
        figureB = figure; set(figureB,'position',imf_pos); %set(figureB,'position',[850 250 800 700]);
        imagesc(IM2.x(1,xmini:xmaxi),IM2.y(1,ymini:ymaxi),IM2.z_masked(ymini:ymaxi,xmini:xmaxi)); axis xy equal; colormap gray; hold on;
        set(gca,'fontsize',14);
        title(['Late date: ',num2str(datestr(date_f,'yyyy/mm/dd'))],'fontsize',16);
    end  
else
    disp('No iceberg coordinate files found. Rerun previous step or check directory.');
end

%plot the DEMs
xmine = nearestneighbour(IM1.x(xminie),DEM1.x); xmaxe = nearestneighbour(IM1.x(xmaxie),DEM1.x);
ymine = nearestneighbour(IM1.y(yminie),DEM1.y); ymaxe = nearestneighbour(IM1.y(ymaxie),DEM1.y);
xmine(xmine<=1) = 1; xmaxe(xmaxe >= length(DEM1.x)) = length(DEM1.x);
ymine(ymine<=1) = 1; ymaxe(ymaxe >= length(DEM1.y)) = length(DEM1.y);
xmin = nearestneighbour(IM2.x(xmini),DEM2.x); xmax = nearestneighbour(IM2.x(xmaxi),DEM2.x);
ymin = nearestneighbour(IM2.y(ymini),DEM2.y); ymax = nearestneighbour(IM2.y(ymaxi),DEM2.y);
xmin(xmin<=1) = 1; xmax(xmax >= length(DEM2.x)) = length(DEM2.x);
ymin(ymin<=1) = 1; ymax(ymax >= length(DEM2.y)) = length(DEM2.y);
%minimum elevations in the plots
zmins = min(DEM1.z_elpsd_adjust(ymine:ymaxe,xmine:xmaxe)); ymins = min(DEM2.z_elpsd_adjust(ymin:ymax,xmin:xmax));
zmins(zmins<-100 | zmins>100) = NaN; ymins(ymins<-100 | ymins>100) = NaN;
zmins(zmins<(nanmedian(zmins)-3*1.4826*mad(zmins,1)) | zmins>(nanmedian(zmins)+3*1.4826*mad(zmins,1))) = NaN;
ymins(ymins<(nanmedian(ymins)-3*1.4826*mad(ymins,1)) | ymins>(nanmedian(ymins)+3*1.4826*mad(ymins,1))) = NaN; 
cmin_early = nanmedian(zmins); cmin_late = nanmedian(ymins);
%set the max for the clim to the 95th percentile of the elevations
figure; h = histogram(DEM1.z_elpsd_adjust(ymine:ymaxe,xmine:xmaxe),min(min(DEM1.z_elpsd_adjust(ymine:ymaxe,xmine:xmaxe))):1:max(max(DEM1.z_elpsd_adjust(ymine:ymaxe,xmine:xmaxe)))); 
zmaxs = h.BinEdges(find(cumsum(h.Values)>=0.9*sum(h.Values),1,'first')+1); clear h; close(gcf);
figure; h = histogram(DEM2.z_elpsd_adjust(ymin:ymax,xmin:xmax),min(min(DEM2.z_elpsd_adjust(ymin:ymax,xmin:xmax))):1:max(max(DEM2.z_elpsd_adjust(ymin:ymax,xmin:xmax)))); 
ymaxs = h.BinEdges(find(cumsum(h.Values)>=0.9*sum(h.Values),1,'first')+1); clear h; close(gcf);
drawnow;

%determine maximum elevation used for color mapping
[ZXgrid,ZYgrid] = meshgrid(DEM1.x,DEM1.y);
xpoly = [early_x-250 early_x+250 early_x+250 early_x-250 early_x-250];
ypoly = [early_y-250 early_y-250 early_y+250 early_y+250 early_y-250];
in = inpolygon(ZXgrid,ZYgrid,xpoly,ypoly);
zmax1 = max(DEM1.z_elpsd_adjust(in))-cmin_early;
clear Z*grid *poly in;
[ZXgrid,ZYgrid] = meshgrid(DEM2.x,DEM2.y);
xpoly = [late_x-250 late_x+250 late_x+250 late_x-250 late_x-250];
ypoly = [late_y-250 late_y-250 late_y+250 late_y+250 late_y-250];
in = inpolygon(ZXgrid,ZYgrid,xpoly,ypoly);
zmax2 = max(DEM2.z_elpsd_adjust(in))-cmin_late;
clear Z*grid *poly in;

%plot DEMs
Znans = isnan(DEM1.z_elpsd_adjust); Ynans = isnan(DEM2.z_elpsd_adjust); 
DEM1.z_elpsd_adjust(Znans) = 0; DEM2.z_elpsd_adjust(Ynans) = 0;
figure1 = figure; set(figure1,'position',DEMo_pos); %set(figure1,'position',[50 800 800 700]);
imagesc(DEM1.x,DEM1.y,DEM1.z_elpsd_adjust-cmin_early); hold on; axis xy equal; 
% set(gca,'clim',[0 max(max(DEM1.z_elpsd_adjust(ymine:ymaxe,xmine:xmaxe)))-cmin_early],'fontsize',14); 
set(gca,'clim',[0 zmax1],'fontsize',14); 
colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
set(gca,'xlim',[min([DEM1.x(xmine);DEM1.x(xmaxe)]) max([DEM1.x(xmine);DEM1.x(xmaxe)])],'ylim',[min([DEM1.y(ymine);DEM1.y(ymaxe)]) max([DEM1.y(ymine);DEM1.y(ymaxe)])]);
plot(early_x,early_y,'*w','markersize',7,'linewidth',1);
title(['Early date: ',num2str(datestr(date_o,'yyyy/mm/dd'))],'fontsize',16); 
figure2 = figure; set(figure2,'position',DEMf_pos); %set(figure2,'position',[50 50 800 700]);
imagesc(DEM2.x,DEM2.y,DEM2.z_elpsd_adjust-cmin_late); hold on; axis xy equal; 
% set(gca,'clim',[0 max(max(DEM2.z_elpsd_adjust(ymin:ymax,xmin:xmax)))-cmin_late],'fontsize',14);
set(gca,'clim',[0 zmax2],'fontsize',14);
colormap(gca,elev_cmap); set(get(cbar,'ylabel'),'string', 'elevation (m)');
set(gca,'xlim',[min([DEM2.x(xmin);DEM2.x(xmax)]) max([DEM2.x(xmin);DEM2.x(xmax)])],'ylim',[min([DEM2.y(ymin);DEM2.y(ymax)]) max([DEM2.y(ymin);DEM2.y(ymax)])]);
plot(early_x,early_y,'*w','markersize',7,'linewidth',1);
title(['Late date: ',num2str(datestr(date_f,'yyyy/mm/dd'))],'fontsize',16); 
figure(figure1);

% make directories for iceberg and icefree ROIs is they do not exist
if ~exist([dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_rois/'],'dir')
    mkdir([dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_rois/']);
end
if ~exist([dir_output,'/',DEM1.time,'-',DEM2.time,'/icefree_rois/'],'dir')
    mkdir([dir_output,'/',DEM1.time,'-',DEM2.time,'/icefree_rois/'])
end

%draw a polygon just inside the iceberg edges in the early DEM
cd([dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_rois/']);
if ~isempty(icebergs_dz)
    for i = 1:length(icebergs_dz); iceberg_list(i,:) = icebergs_dz(i).name(8:9); end
else
    iceberg_list = 'XX';
end
if isempty(strmatch(iceberg_no,iceberg_list))
    %draw the iceberg ROI
    figure(figure1);
    disp('Zoom in on the iceberg (leaving a buffer containing thin ice or open water) in the early DEM by clicking the UL & LR corners');
    [a] = ginput(2); %get the UL & LR corner coordinates
    xref = nearestneighbour(a(:,1)',DEM1.x); yref = nearestneighbour(a(:,2)',DEM1.y);
    set(gca,'xlim',[min(DEM1.x(xref)) max(DEM1.x(xref))],'ylim',[min(DEM1.y(yref)) max(DEM1.y(yref))]);
    clear xref yref;
    figure(figureA);
    xref = nearestneighbour(a(:,1)',IM1.x); yref = nearestneighbour(a(:,2)',IM1.y);
%     set(gca,'xlim',[min(im1.x(xref)) max(im1.x(xref))],'ylim',[min(im1.y(yref)) max(im1.y(yref))]);
    set(gca,'xlim',[min(IM1.x(xref)) max(IM1.x(xref))],'ylim',[min(IM1.y(yref)) max(IM1.y(yref))]);
    clear xref yref;
    figure(figure1);
    disp('When crosshairs appear on the early DEM, click on vertices just inside the iceberg edge to draw a polygon');
    disp('NOTE: white areas = DEM holes');
    [~,xm,ym] = roipoly;
    xmi = nearestneighbour(xm',DEM1.x); ymi = nearestneighbour(ym',DEM1.y);
    plot(DEM1.x(xmi),DEM1.y(ymi),'-*k','linewidth',2,'markersize',4); hold on;
%     set(gca,'xlim',[min(DEM1.x(xmi))-DEMhalo1 max(DEM1.x(xmi))+DEMhalo1],'ylim',[min(DEM1.y(ymi))-DEMhalo1 max(DEM1.y(ymi))+DEMhalo1]);
    S.Geometry = 'Polygon';
    S.BoundingBox = [min(xm) min(ym); max(xm) max(ym)];
    S.X = double(xm'); S.Y = double(ym');
    S.Name = ['iceberg',iceberg_no];
    shapefile_name = [region_abbrev,'_',num2str(DEM1.time),'_iceberg',iceberg_no];
    shapewrite(S,shapefile_name);
    copyfile([dir_code,PSprojfile],[dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_rois/',shapefile_name,'.prj']);
    clear xm ym;
    figure(figureA);
%     set(gca,'xlim',[min([DEM1.x(xmine);DEM1.x(xmaxe)]) max([DEM1.x(xmine);DEM1.x(xmaxe)])],'ylim',[min([DEM1.y(ymine);DEM1.y(ymaxe)]) max([DEM1.y(ymine);DEM1.y(ymaxe)])]);
    
    %draw the icefree ROI in the early image
    cd([dir_output,'/',DEM1.time,'-',DEM2.time,'/icefree_rois/']);
    figure(figure1);
    disp('Now click on vertices of a nearby low elevation region to draw an initial guess for an ice-free polygon');
    [~,xm,ym] = roipoly;
    xmif = nearestneighbour(xm',DEM1.x); ymif = nearestneighbour(ym',DEM1.y);
    plot(DEM1.x(xmif),DEM1.y(ymif),'-*y','linewidth',2,'markersize',4); hold on;
    SL.Geometry = 'Polygon';
    SL.BoundingBox = [min(xm) min(ym); max(xm) max(ym)];
    SL.X = double(xm'); SL.Y = double(ym');
    SL.Name = ['icefree',iceberg_no];
    shapefile_name = [region_abbrev,'_',num2str(DEM1.time),'_icefree',iceberg_no];
    shapewrite(SL,shapefile_name);
    copyfile([dir_code,PSprojfile],[dir_output,'/',DEM1.time,'-',DEM2.time,'/icefree_rois/',shapefile_name,'.prj']);
    sl_early.X = SL.X(1:size(SL.X,2)); sl_early.Y = SL.Y(1:size(SL.Y,2)); %truncate last value (NaN)
    clear SL xm ym xmif ymif;
    
    %assign a quality flag to the ice-free polygon
    prompt = 'Is open water or thin ice present (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        sl_early.quality = 1; %higher quality adjustments assigned a value = 1
    else
        sl_early.quality = 0; %lower quality adjustments assigned a value = 0
    end
    
    %draw the icefree ROI in the later image
    figure(figure2);
    disp('Now zoom in on the iceberg (leaving a buffer containing thin ice or open water) in the later DEM by clicking the UL & LR corners');
    [a] = ginput(2); %get the UL & LR corner coordinates
    xref = nearestneighbour(a(:,1)',DEM2.x); yref = nearestneighbour(a(:,2)',DEM2.y);
    set(gca,'xlim',[min(DEM2.x(xref)) max(DEM2.x(xref))],'ylim',[min(DEM2.y(yref)) max(DEM2.y(yref))]);
    clear xref yref;
    figure(figureB);
    xref = nearestneighbour(a(:,1)',IM2.x); yref = nearestneighbour(a(:,2)',IM2.y);
%     set(gca,'xlim',[min(im1.x(xref)) max(im1.x(xref))],'ylim',[min(im1.y(yref)) max(im1.y(yref))]);
    set(gca,'xlim',[min(IM2.x(xref)) max(IM2.x(xref))],'ylim',[min(IM2.y(yref)) max(IM2.y(yref))]);
    clear xref yref;
%     figure(figureB); set(gca,'xlim',[min([DEM2.x(xmin);DEM2.x(xmax)]) max([DEM2.x(xmin);DEM2.x(xmax)])],'ylim',[min([DEM2.y(ymin);DEM2.y(ymax)]) max([DEM2.y(ymin);DEM2.y(ymax)])]);
    figure(figure2);
    disp('When crosshairs appear on the later DEM, click on vertices of a nearby low elevation region to draw an initial guess for an ice-free polygon');
    [~,xm,ym] = roipoly;
    xmif = nearestneighbour(xm',DEM2.x); ymif = nearestneighbour(ym',DEM2.y);
    plot(DEM2.x(xmif),DEM2.y(ymif),'-*y','linewidth',2,'markersize',4); hold on;
    SL.Geometry = 'Polygon';
    SL.BoundingBox = [min(xm) min(ym); max(xm) max(ym)];
    SL.X = double(xm'); SL.Y = double(ym');
    SL.Name = ['icefree',iceberg_no];
    shapefile_name = [region_abbrev,'_',num2str(DEM2.time),'_icefree',iceberg_no];
    shapewrite(SL,shapefile_name);
    copyfile([dir_code,PSprojfile],[dir_output,'/',DEM1.time,'-',DEM2.time,'/icefree_rois/',shapefile_name,'.prj']);
    sl_late.X = SL.X(1:size(SL.X,2)); sl_late.Y = SL.Y(1:size(SL.Y,2)); %truncate last value (NaN)
    clear SL xm ym xmif ymif;
    
    %assign a quality flag to the ice-free polygon
    prompt = 'Is open water or thin ice present (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        sl_late.quality = 1; %higher quality adjustments assigned a value = 1
    else
        sl_late.quality = 0; %lower quality adjustments assigned a value = 0
    end
    
else
    disp('Use existing shapefiles');
    iceberg_file = [dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_rois/',region_abbrev,'_',num2str(DEM1.time),'_iceberg',iceberg_no,'.shp'];
    S = shaperead(iceberg_file);
    S.X = S.X(1:size(S.X,2)-1); S.Y = S.Y(1:size(S.Y,2)-1); %truncate last value (NaN)
    
    %plot the iceberg polygon on the earlier DEM
    figure(figure1);
    plot(S.X,S.Y,'-*k','linewidth',2,'markersize',4); hold on;
    drawnow;
    
    %read the shape-file for nearby ice-free pixels in the early DEM
    icefree_early_file = [dir_output,'/',DEM1.time,'-',DEM2.time,'/icefree_rois/',region_abbrev,'_',num2str(DEM1.time),'_icefree',iceberg_no,'.shp'];
    cd_to_icefreerois = ['cd ',dir_output,'/',DEM1.time,'-',DEM2.time,'/icefree_rois/']; eval(cd_to_icefreerois);
    sl_early = shaperead(icefree_early_file); 
    sl_early.X = sl_early.X(1:size(sl_early.X,2)-1); sl_early.Y = sl_early.Y(1:size(sl_early.Y,2)-1); %truncate last value (NaN)
    figure(figureA); set(gca,'xlim',[min([DEM1.x(xmine);DEM1.x(xmaxe)]) max([DEM1.x(xmine);DEM1.x(xmaxe)])],'ylim',[min([DEM1.y(ymine);DEM1.y(ymaxe)]) max([DEM1.y(ymine);DEM1.y(ymaxe)])]);
    prompt = 'Is open water or thin ice present in the early image (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        sl_early.quality = 1; %higher quality adjustments assigned a value = 1
    else
        sl_early.quality = 0; %lower quality adjustments assigned a value = 0
    end
    icefree_late_file = [region_abbrev,'_',num2str(DEM2.time),'_icefree',iceberg_no,'.shp'];
    sl_late = shaperead(icefree_late_file); cd ../iceberg_rois
    sl_late.X = sl_late.X(1:size(sl_late.X,2)-1); sl_late.Y = sl_late.Y(1:size(sl_late.Y,2)-1); %truncate last value (NaN)
    figure(figureB); set(gca,'xlim',[min([DEM2.x(xmin);DEM2.x(xmax)]) max([DEM2.x(xmin);DEM2.x(xmax)])],'ylim',[min([DEM2.y(ymin);DEM2.y(ymax)]) max([DEM2.y(ymin);DEM2.y(ymax)])]);
    prompt = 'Is open water or thin ice present in the late image (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        sl_late.quality = 1; %higher quality adjustments assigned a value = 1
    else
        sl_late.quality = 0; %lower quality adjustments assigned a value = 0
    end
    
    %if necessary, re-draw the iceberg polygon in the early DEM
    figure(figure1);
    disp('click on the early DEM to activate');
    figure(figure1); ax = waitforbuttonpress; clear ax;
    disp('Zoom in on the iceberg in the DEM by clicking the UL & LR corners');
    [a] = ginput(2); %get the UL & LR corner coordinates
    xref = nearestneighbour(a(:,1)',DEM1.x); yref = nearestneighbour(a(:,2)',DEM1.y);
    set(gca,'xlim',[min(DEM1.x(xref)) max(DEM1.x(xref))],'ylim',[min(DEM1.y(yref)) max(DEM1.y(yref))]);
    clear xref yref;
    figure(figureA);
%         xref = nearestneighbour(a(:,1)',im1.x); yref = nearestneighbour(a(:,2)',im1.y);
%         set(gca,'xlim',[min(im1.x(xref)) max(im1.x(xref))],'ylim',[min(im1.y(yref)) max(im1.y(yref))]);
    xref = nearestneighbour(a(:,1)',IM1.x); yref = nearestneighbour(a(:,2)',IM1.y);
    set(gca,'xlim',[min(IM1.x(xref)) max(IM1.x(xref))],'ylim',[min(IM1.y(yref)) max(IM1.y(yref))]);
    clear xref yref;

    prompt = 'Look at the mask overlay on the DEM. Does it need to be modified (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        clear S;
        figure(figure1);
        disp('click on the early DEM to activate');
        figure(figure1); ax = waitforbuttonpress; clear ax;
        disp('When crosshairs appear, click on iceberg vertices to draw a polygon');
        [~,xm,ym] = roipoly;
        xmi = nearestneighbour(xm',DEM1.x); ymi = nearestneighbour(ym',DEM1.y);
        plot(DEM1.x(xmi),DEM1.y(ymi),'-*k','linewidth',2,'markersize',4); hold on;

        cd([dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_rois/']);
        S.Geometry = 'Polygon';
        S.BoundingBox = [min(xm) min(ym); max(xm) max(ym)];
        S.X = double(xm'); S.Y = double(ym');
        S.Name = ['iceberg',iceberg_no];
        shapefile_name = [region_abbrev,'_',num2str(DEM1.time),'_iceberg',iceberg_no];
        copyfile([dir_code,PSprojfile],[dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_rois/',shapefile_name,'.prj']);
        cd([dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_rois/']);
    else
        xmi = nearestneighbour(S.X,DEM1.x); ymi = nearestneighbour(S.Y,DEM1.y);
    end
end

%use the image to tweek the iceberg mask if necessary
figure(figure1); figure(figureA);
plot(DEM1.x(xmi),DEM1.y(ymi),'-*r','linewidth',2,'markersize',4); hold on;
prompt = 'Based on the mask overlay on the early IMAGE, does the iceberg mask need to be modified (y/n)?';
str = input(prompt,'s');
if strmatch(str,'y')==1
    clear S; figure(figureA); set(figureA,'position',imo_pos); %set(figureA,'position',[450 50 1200 1100]);
    clear xmi ymi;
    disp('click on the early DEM to activate');
    figure(figure1); ax = waitforbuttonpress; clear ax;
    disp('Zoom in on the iceberg in the early DEM by clicking the UL & LR corners');
    [a] = ginput(2); %get the UL & LR corner coordinates
    xref = nearestneighbour(a(:,1)',DEM1.x); yref = nearestneighbour(a(:,2)',DEM1.y);
    set(gca,'xlim',[min(DEM1.x(xref)) max(DEM1.x(xref))],'ylim',[min(DEM1.y(yref)) max(DEM1.y(yref))]);
    clear xref yref;
    figure(figureA);
    xref = nearestneighbour(a(:,1)',IM1.x); yref = nearestneighbour(a(:,2)',IM1.y);
%     set(gca,'xlim',[min(im1.x(xref)) max(im1.x(xref))],'ylim',[min(im1.y(yref)) max(im1.y(yref))]);
    set(gca,'xlim',[min(IM1.x(xref)) max(IM1.x(xref))],'ylim',[min(IM1.y(yref)) max(IM1.y(yref))]);
    clear xref yref;
%     clear xref yref; clear a;
    disp('When crosshairs appear, click on iceberg vertices to draw a polygon');
    [~,xm,ym] = roipoly;
    xmi = nearestneighbour(xm',DEM1.x); ymi = nearestneighbour(ym',DEM1.y);
    plot(DEM1.x(xmi),DEM1.y(ymi),'--*','color',[0.5 0.5 0.5],'linewidth',2,'markersize',4); hold on;
    set(figureA,'position',imo_pos); %set(figureA,'position',[850 800 800 700]);
    set(gca,'xlim',[min([DEM1.x(xmine);DEM1.x(xmaxe)]) max([DEM1.x(xmine);DEM1.x(xmaxe)])],'ylim',[min([DEM1.y(ymine);DEM1.y(ymaxe)]) max([DEM1.y(ymine);DEM1.y(ymaxe)])]);
    figure(figure1);
    plot(DEM1.x(xmi),DEM1.y(ymi),'--*','color',[0.5 0.5 0.5],'linewidth',2,'markersize',4); hold on;
    set(gca,'xlim',[min([DEM1.x(xmine);DEM1.x(xmaxe)]) max([DEM1.x(xmine);DEM1.x(xmaxe)])],'ylim',[min([DEM1.y(ymine);DEM1.y(ymaxe)]) max([DEM1.y(ymine);DEM1.y(ymaxe)])]);
    drawnow;
    
    cd([dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_rois/']);
    S.Geometry = 'Polygon';
    S.BoundingBox = [min(xm) min(ym); max(xm) max(ym)];
    S.X = double(xm'); S.Y = double(ym');
    S.Name = ['iceberg',iceberg_no];
    shapefile_name = [region_abbrev,'_',num2str(DEM1.time),'_iceberg',iceberg_no];
    shapewrite(S,shapefile_name);
    copyfile([dir_code,PSprojfile],[dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_rois/',shapefile_name,'.prj']);
    cd([dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_rois/']);
end
drawnow;
disp('Done modifying iceberg polygon for elevation extraction');

%calculate the iceberg dimensions to crop the DEMs: use the max of both
%dimensions to account for possible rotation of big skinny icebergs
berg_width = abs(max(xmi) - min(xmi));
berg_length = abs(max(ymi) - min(ymi));
xmine = xce-(max([berg_width,berg_length])+50); xmaxe = xce+(max([berg_width,berg_length])+50);
ymine = yce-(max([berg_width,berg_length])+50); ymaxe = yce+(max([berg_width,berg_length])+50);
xmine(xmine<=1) = 1; xmaxe(xmaxe >= length(DEM1.x)) = length(DEM1.x);
ymine(ymine<=1) = 1; ymaxe(ymaxe >= length(DEM1.y)) = length(DEM1.y);
xmin = xc-(max([berg_width,berg_length])+50); xmax = xc+(max([berg_width,berg_length])+50);
ymin = yc-(max([berg_width,berg_length])+50); ymax = yc+(max([berg_width,berg_length])+50);
xmin(xmin<=1) = 1; xmax(xmax >= length(DEM2.x)) = length(DEM2.x);
ymin(ymin<=1) = 1; ymax(ymax >= length(DEM2.y)) = length(DEM2.y);

%crop dataset extents
disp('Cropping DEMs to speed-up processing time');
A.x = DEM1.x(1,xmine:xmaxe); A.y = DEM1.y(1,ymine:ymaxe);
B.x = DEM2.x(1,xmin:xmax); B.y = DEM2.y(1,ymin:ymax);
A.z = DEM1.z(ymine:ymaxe,xmine:xmaxe);
B.z = DEM2.z(ymin:ymax,xmin:xmax);
early_elev_offset = DEM1.geoid_z(ymine:ymaxe,xmine:xmaxe);
late_elev_offset = DEM2.geoid_z(ymin:ymax,xmin:xmax);
DEM1.z_elpsd_adjust(Znans) = NaN; DEM2.z_elpsd_adjust(Ynans) = NaN;
A.z_elpsd_adjust = DEM1.z_elpsd_adjust(ymine:ymaxe,xmine:xmaxe); 
% A.z_elpsd_adjust(A.z_elpsd_adjust<0) = 0;
B.z_elpsd_adjust = DEM2.z_elpsd_adjust(ymin:ymax,xmin:xmax); 
% B.z_elpsd_adjust(B.z_elpsd_adjust<0) = 0;

%clear DEM1 & DEM2 structures to speed-up computation time
%clear DEM1 DEM2;

%estimate fjord offset from regional coregistration
fjordz_o = fjord.dem1.z(~isnan(fjord.dem1.z)); fjordz_o(fjordz_o>100) = NaN; fjordz_o(fjordz_o<-100) = NaN; 
fjordz_medo = median(fjord.dem1.z(~isnan(fjord.dem1.z)),'omitnan'); fjordz_mado = mad(fjord.dem1.z(~isnan(fjord.dem1.z)),1);
fjordz_o(fjordz_o<fjordz_medo-3*(1.4826*fjordz_mado) | fjordz_o>fjordz_medo+3*(1.4826*fjordz_mado)) = NaN;
fjordz_f = fjord.dem2.z(~isnan(fjord.dem2.z)); fjordz_f(fjordz_f>100) = NaN; fjordz_f(fjordz_f<-100) = NaN; 
fjordz_medf = median(fjord.dem2.z(~isnan(fjord.dem2.z)),'omitnan'); fjordz_madf = mad(fjord.dem2.z(~isnan(fjord.dem2.z)),1);
fjordz_f(fjordz_f<fjordz_medf-3*(1.4826*fjordz_madf) | fjordz_f>fjordz_medf+3*(1.4826*fjordz_madf)) = NaN;
fjord_offset = nanmedian(min(fjordz_o)) - nanmedian(min(fjordz_f));

%(not recommended) uncomment below to enable use of bedrock & tide-modelled offset for regional sea level
%adjustment if tidal model was used & this is desired
% dtide = fjord.tidal_change;
% bed_offset_file = dir([dir_output,'bedrock_offset*.mat']);
% if ~isempty(bed_offset_file) %if overlapping bedrock regions were used to adjust DEMs
%     load([dir_output,'bedrock_offset_',num2str(DEM1.time),'-',num2str(DEM2.time),'.mat']);
%     
%     %calculate the median bedrock offset between DEMs
%     reg1_dz = []; reg2_dz = [];
%     for i = 1:size(bedrock.reg1.dem_diff.map,1)
%         for j = 1:size(bedrock.reg1.dem_diff.map,2)
%             reg1_dz(size(reg1_dz,2)+1) = bedrock.reg1.dem_diff.map(i,j);
%         end
%     end
%     for i = 1:size(bedrock.reg2.dem_diff.map,1)
%         for j = 1:size(bedrock.reg2.dem_diff.map,2)
%             reg2_dz(size(reg2_dz,2)+1) = bedrock.reg2.dem_diff.map(i,j);
%         end
%     end
%     bias = nanmedian([reg1_dz(~isnan(reg1_dz)) reg2_dz(~isnan(reg2_dz))]);
%     fjord_offset(1) = nanmedian(min(fjordz_o)) - nanmedian(min(fjordz_f));
%     fjord_offset(2) = bias + dtide;
% end
A.z_reg_adjust = A.z_elpsd_adjust - nanmedian(min(fjordz_o)); B.z_reg_adjust = B.z_elpsd_adjust - nanmedian(min(fjordz_f));

%pull max elevation estimate from the iceberg
poly1 = roipoly(A.x,A.y,A.z_reg_adjust,S.X,S.Y);
mask1 = logical(poly1);
berg_elevs = mask1.*A.z_reg_adjust; berg_elevs(berg_elevs == 0) = NaN;
% cmax=ceil(max(max(berg_elevs)));
figure; h = histogram(berg_elevs,min(berg_elevs(~isnan(berg_elevs))):1:max(berg_elevs(~isnan(berg_elevs)))); 
cmax = h.BinEdges(find(cumsum(h.Values)>=0.98*sum(h.Values),1,'first')+1); clear h; close(gcf); drawnow;
figure(figure1); set(gca,'clim',[0 cmax]);
figure(figure2); set(gca,'clim',[0 cmax]);

%remove DEM bias using local ice-free regions
prompt = 'Are there good water/sea ice elevations around the iceberg in both DEMs (y/n)?'; ibmask_str = input(prompt,'s');
if strmatch(ibmask_str,'y')==1
    [A,B,sl_offset1,sl_offset2,cmin_early,cmin_late,flag] = sea_level_adjust(DEM1.time,DEM2.time,IM1,IM2,dir_output,A,B,S,sl_early,sl_late,fjord_offset,fjordz_o,fjordz_f,iceberg_no,region_abbrev,cmax,DEMo_pos,DEMf_pos,imo_pos,imf_pos);
else
    disp('using regional sea level adjustment for elevation change estimates');
    sl_offset1 = nanmedian(min(fjordz_o)); sl_offset2 = nanmedian(min(fjordz_f));
    flag = 1; %note that the adjustment is not reliable
    
    %replot figures with regional elevation adjustments for coregistration
    clf(figure1); figure(figure1);
    A.z_local_adjust = A.z_elpsd_adjust - sl_offset1;
    imagesc(A.x,A.y,A.z_reg_adjust); hold on; axis xy equal; set(gca,'clim',[0 (cmax-(sl_offset2 + sl_offset1)/2)],'fontsize',14);
    cmap = colormap(gca,elev_cmap); cmap(1,:) = [1 1 1]; colormap(gca,cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
    title(['Early date: ',num2str(DEM1.time)],'fontsize',16);
    plot(S.X,S.Y,'-*k','linewidth',2,'markersize',4); hold on; %iceberg ROI
    B.z_local_adjust = B.z_elpsd_adjust - sl_offset2;
    clf(figure2); figure(figure2);
    imagesc(B.x,B.y,B.z_reg_adjust); hold on; axis xy equal; set(gca,'clim',[0 (cmax-(sl_offset2 + sl_offset1)/2)],'fontsize',14);
    cmap = colormap(gca,elev_cmap); cmap(1,:) = [1 1 1]; colormap(gca,cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
    title(['Late date: ',num2str(DEM2.time)],'fontsize',16);
    
    zmins2 = min(A.z_local_adjust); ymins2 = min(B.z_local_adjust);
    zmins2(zmins2<-100 | zmins2>100) = NaN; ymins2(ymins2<-100 | ymins2>100) = NaN;
    zmins2(zmins2<(nanmedian(zmins2)-3*1.4826*mad(zmins2,1)) | zmins2>(nanmedian(zmins2)+3*1.4826*mad(zmins2,1))) = NaN;
    ymins2(ymins2<(nanmedian(ymins2)-3*1.4826*mad(ymins2,1)) | ymins2>(nanmedian(ymins2)+3*1.4826*mad(ymins2,1))) = NaN;
    cmin_early = floor(nanmedian(zmins2)); cmin_late = floor(nanmedian(ymins2));
end
cmin = min([cmin_early;cmin_late]);

%pull updated max elevation estimate from the iceberg in the local elevation-adjusted DEM
poly1 = roipoly(A.x,A.y,A.z_local_adjust,S.X,S.Y);
mask1 = logical(poly1);
berg_elevs = mask1.*A.z_local_adjust; berg_elevs(berg_elevs == 0) = NaN;
cmax=ceil(max(max(berg_elevs)));

%crop the area around the iceberg of interest in the early DEM
disp('Zoom in on iceberg & crop dataset in earlier DEM by clicking on its upper left & lower right corners (leave a buffer containing ice-free data!)');
close all; drawnow;
Anans = isnan(A.z_local_adjust); Bnans = isnan(B.z_local_adjust); 
A.z_local_adjust(Anans) = -50; B.z_local_adjust(Bnans) = -50;
figure1 = figure; set(figure1,'position',DEMo_pos); %set(figure1,'position',[50 250 800 700]); 
clear a;
imagesc(A.x,A.y,A.z_local_adjust); hold on; axis xy equal; set(gca,'clim',[0 cmax],'fontsize',14);
colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
title(['Early date: ',num2str(DEM1.time)],'fontsize',16);
plot(S.X,S.Y,'-*k','linewidth',2,'markersize',4); hold on;
[a] = ginput(2); %get the UL & LR corner coordinates
xref = nearestneighbour(a(:,1)',A.x); yref = nearestneighbour(a(:,2)',A.y);
C.x = A.x(1,min(xref):max(xref)); C.y = A.y(1,min(yref):max(yref));
C.z = A.z(min(yref):max(yref),min(xref):max(xref));
C.z_elpsd_adjust = A.z_elpsd_adjust(min(yref):max(yref),min(xref):max(xref));
C.z_local_adjust = A.z_local_adjust(min(yref):max(yref),min(xref):max(xref));
C.geoid_z = early_elev_offset(min(yref):max(yref),min(xref):max(xref));
clear A; A = C; clear C; 
clf(figure1); figure(figure1);
imagesc(A.x,A.y,A.z_local_adjust); hold on; axis xy equal; set(gca,'clim',[cmin cmax],'fontsize',14);
colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
title(['Early date: ',num2str(DEM1.time)],'fontsize',16);
plot(S.X,S.Y,'-*k','linewidth',2,'markersize',4); hold on; drawnow;
drawnow;

%crop the area around the iceberg of interest in the later DEM
disp('Zoom in on iceberg & crop dataset in later DEM by clicking on its upper left & lower right corners (leave a buffer containing ice-free data!)');
figure2 = figure; set(figure2,'position',DEMf_pos); %set(figure2,'position',[850 250 800 700]); 
clear a xref yref;
imagesc(B.x,B.y,B.z_local_adjust); hold on; axis xy equal; set(gca,'clim',[cmin cmax],'fontsize',14);
colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
title(['Late date: ',num2str(DEM2.time)],'fontsize',16);
[a] = ginput(2);
xref = nearestneighbour(a(:,1)',B.x); yref = nearestneighbour(a(:,2)',B.y);
C.x = B.x(1,min(xref):max(xref)); C.y = B.y(1,min(yref):max(yref));
C.z = B.z(min(yref):max(yref),min(xref):max(xref));
C.z_elpsd_adjust = B.z_elpsd_adjust(min(yref):max(yref),min(xref):max(xref));
C.z_local_adjust = B.z_local_adjust(min(yref):max(yref),min(xref):max(xref));
C.geoid_z = late_elev_offset(min(yref):max(yref),min(xref):max(xref));
clear B; B = C; clear C; 
clf(figure2); figure(figure2);
imagesc(B.x,B.y,B.z_local_adjust); hold on; axis xy equal; set(gca,'clim',[0 cmax],'fontsize',14);
colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
title(['Late date: ',num2str(DEM2.time)],'fontsize',16);
drawnow;

%plot contours over the iceberg surface to help with feature ID
figure(figure1);
[cont,conth] = contour(A.x,A.y,A.z_local_adjust,[0:1:round(cmax)]);
conth.LineColor = 'k';
figure(figure2);
[cont,conth] = contour(B.x,B.y,B.z_local_adjust,[0:1:round(cmax)]);
conth.LineColor = 'k';

%mask-out regions within icebergs that have anomalous elevations (holes or spikes)
disp('If spurious elevations (big spikes or holes) are evident in the iceberg DEMs, draw masks around them to remove them');
p=1;
while p
    if p == 1
        prompt = 'Create mask for early DEM (y/n)?'; str = input(prompt,'s');
        spurious_o = zeros(size(A.z_local_adjust));
    else
        prompt = 'Are additional masks necessary (y/n)?'; str = input(prompt,'s');
    end
    if strmatch(str,'y')==1
        figure(figure1);
        anom_z_mask1 = roipoly; anom_z_mask1 = double(~anom_z_mask1);
        A.z = anom_z_mask1.*A.z; A.z(anom_z_mask1==0) = NaN;
        A.z_elpsd_adjust = anom_z_mask1.*A.z_elpsd_adjust; A.z_elpsd_adjust(anom_z_mask1==0) = NaN;
        A.z_local_adjust = anom_z_mask1.*A.z_local_adjust; A.z_local_adjust(anom_z_mask1==0) = NaN;
        imagesc(A.x,A.y,A.z_local_adjust); hold on; axis xy equal; set(gca,'clim',[cmin cmax],'fontsize',14);
        colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
        title(['Early date: ',num2str(DEM1.time)],'fontsize',16);
        plot(S.X,S.Y,'-*k','linewidth',2,'markersize',4); hold on; 
        [cont,conth] = contour(A.x,A.y,A.z_local_adjust,[0:1:round(cmax)]);
        conth.LineColor = 'k';
        drawnow;
        blunders = ~anom_z_mask1;
        spurious_o = spurious_o + blunders;
        p = p+1;
    else
        spurious_o(spurious_o>0) = 1;
        break
    end
end
p=1;
while p
    if p == 1
        prompt = 'Create mask for later DEM (y/n)?'; str = input(prompt,'s');
        spurious_f = zeros(size(B.z_local_adjust));
    else
        prompt = 'Are additional masks necessary (y/n)?'; str = input(prompt,'s');
    end
    if strmatch(str,'y')==1
        figure(figure2);
        anom_z_mask2 = roipoly; anom_z_mask2 = double(~anom_z_mask2);
        B.z = anom_z_mask2.*B.z; B.z(anom_z_mask2==0) = NaN;
        B.z_elpsd_adjust = anom_z_mask2.*B.z_elpsd_adjust; B.z_elpsd_adjust(anom_z_mask2==0) = NaN;
        B.z_local_adjust = anom_z_mask2.*B.z_local_adjust; B.z_local_adjust(anom_z_mask2==0) = NaN;
        imagesc(B.x,B.y,B.z_local_adjust); hold on; hold on; axis xy equal; set(gca,'clim',[cmin cmax],'fontsize',14);
        colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
        title(['Late date: ',num2str(DEM2.time)],'fontsize',16);
        [cont,conth] = contour(B.x,B.y,B.z_local_adjust,[0:1:round(cmax)]);
        conth.LineColor = 'k';
        drawnow;
        blunders = ~anom_z_mask2;
        spurious_f = spurious_f + blunders;
        p = p+1;
    else
        spurious_f(spurious_f>0) = 1;
        break
    end
end
w=who;%create a variable that lists all variables in the current workspace

%create masks for data gaps
gapo_mask = A.z_elpsd_adjust; gapo_mask(gapo_mask~=0) = 1; gapo_mask(isnan(A.z_elpsd_adjust)) = 0;
gapf_mask = B.z_elpsd_adjust; gapf_mask(gapf_mask~=0) = 1; gapf_mask(isnan(B.z_elpsd_adjust)) = 0;

%recreate elevation maps with the regional adjustment (somehow getting removed)
A.z_reg_adjust = A.z_elpsd_adjust - nanmedian(min(fjordz_o)); B.z_reg_adjust = B.z_elpsd_adjust - nanmedian(min(fjordz_f));

% perform iceberg translation & rotation 5 times to quantify user uncertainty
disp('Extract iceberg elevations from the same portion of the iceberg on each date. Repeat same procedure 5 times!');
p=1;
while p
    disp(['Loop #',num2str(p)]);
    if p == 1
        close all; drawnow;
    end
    
    %determine iceberg translation and rotation between DEM acquisition dates
    [b,c,vertex_dist,vertex_ang,berg_dist,berg_xoffset,berg_yoffset,k] = measure_iceberg_motion(DEM1.time,DEM2.time,IM1,IM2,A,B,S,cmin,cmax,DEMo_pos,DEMf_pos,imo_pos,imf_pos);
    
    %save new vertex coordinates
    disp('Extracting elevation change information from within the iceberg polygon');
    IB(p).to = date_o; IB(p).tf = date_f;
    IB(p).elpsd_adjust_o = -early_elev_offset; IB(p).elpsd_adjust_f = -late_elev_offset;
    IB(p).local_adjust_o = -sl_offset1; IB(p).local_adjust_f = -sl_offset2;
    IB(p).waterflag_o = sl_early.quality; IB(p).waterflag_f = sl_late.quality; 
    IB(p).flag = flag;
    IB(p).vertices.xo = S.X; IB(p).vertices.yo = S.Y;
    IB(p).vertices.xo_mid = b(1,1); IB(p).vertices.yo_mid = b(1,2);
    IB(p).vertices.ango = vertex_ang;
    IB(p).travel_xy_dist = berg_dist;
    IB(p).travel_x_dist = berg_xoffset; IB(p).travel_y_dist = berg_yoffset;
    IB(p).rotate = -k; IB(p).vertices.angf = vertex_ang + IB(p).rotate;
    IB(p).vertices.xf = vertex_dist.*cosd(IB(p).vertices.angf)+c(1,1);
    IB(p).vertices.yf = vertex_dist.*sind(IB(p).vertices.angf)+c(1,2);
    IB(p).vertices.xf_mid = c(1,1); IB(p).vertices.yf_mid = c(1,2);
    IB(p).zo.blunder_mask = spurious_o; IB(p).zf.blunder_mask = spurious_f;
    IB(p).zo.gap_mask = gapo_mask; IB(p).zf.gap_mask = gapf_mask;
    
    %extract elevation change estimates
    clear Anans Bnans;
    Anans = isnan(A.z_elpsd_adjust); Bnans = isnan(B.z_elpsd_adjust); 
    A.z_reg_adjust(Anans) = NaN; B.z_reg_adjust(Bnans) = NaN;
    A.z_local_adjust(Anans) = NaN; B.z_local_adjust(Bnans) = NaN;
    if (p == 1 && ~isempty(strmatch('Sl',w))) || isempty(strmatch('Sl',w))
        disp('Rotating & shifting the iceberg from its position in the later date to the early date location & orientation');
        
        %replace all values outside the iceberg polygon with NaNs
        poly1 = roipoly(A.x,A.y,A.z_local_adjust,IB(p).vertices.xo,IB(p).vertices.yo);
        mask1 = logical(poly1);
        zo.reg_adjust.map = mask1.*A.z_reg_adjust; zo.reg_adjust.map(zo.reg_adjust.map == 0) = NaN; zo.reg_adjust.map(IB(p).zo.gap_mask==0) = NaN;
        zo.local_adjust.map = mask1.*A.z_local_adjust; zo.local_adjust.map(zo.local_adjust.map == 0) = NaN; zo.local_adjust.map(IB(p).zo.gap_mask==0) = NaN;
        poly2 = roipoly(B.x,B.y,B.z_local_adjust,IB(p).vertices.xf,IB(p).vertices.yf);
        mask2 = logical(poly2);
        zf.reg_adjust.map = mask2.*B.z_reg_adjust; zf.reg_adjust.map(zf.reg_adjust.map == 0) = NaN; zf.reg_adjust.map(IB(p).zf.gap_mask==0) = NaN;
        zf.local_adjust.map = mask2.*B.z_local_adjust; zf.local_adjust.map(zf.local_adjust.map == 0) = NaN; zf.local_adjust.map(IB(p).zf.gap_mask==0) = NaN;
        
        %crop the x,y, & masked-z extent to only the iceberg of interest & immediate surroundings
        xo_index = nearestneighbour(IB(p).vertices.xo,A.x); yo_index = nearestneighbour(IB(p).vertices.yo,A.y);
        IB(p).xo = A.x(1,min(xo_index):max(xo_index)); IB(p).yo = A.y(1,min(yo_index):max(yo_index));
        IB(p).zo.reg_adjust.map = zo.reg_adjust.map(min(yo_index):max(yo_index),min(xo_index):max(xo_index));
        IB(p).zo.local_adjust.map = zo.local_adjust.map(min(yo_index):max(yo_index),min(xo_index):max(xo_index));
        xf_index = nearestneighbour(IB(p).vertices.xf,B.x); yf_index = nearestneighbour(IB(p).vertices.yf,B.y);
        IB(p).xf = B.x(1,min(xf_index):max(xf_index)); IB(p).yf = B.y(1,min(yf_index):max(yf_index));
        IB(p).zf.reg_adjust.map = zf.reg_adjust.map(min(yf_index):max(yf_index),min(xf_index):max(xf_index));
        IB(p).zf.local_adjust.map = zf.local_adjust.map(min(yf_index):max(yf_index),min(xf_index):max(xf_index));
        disp('Masked iceberg in (regional and local) sea level-adjusted DEMs');
        
        %translate & rotate the iceberg in the later DEM onto the same earlier reference frame
        [x_rot,y_rot] = unrotate_untranslate_iceberg(IB,p);
        
        %find the nearest neighboring elevation from the earlier DEM to difference
        %with the rotated, translated later DEM
        [xgrid,ygrid] = meshgrid(IB(p).xo,IB(p).yo);
        zfreg_rot = griddata(x_rot,y_rot,double(IB(p).zf.reg_adjust.map),xgrid,ygrid);
        zfloc_rot = griddata(x_rot,y_rot,double(IB(p).zf.local_adjust.map),xgrid,ygrid);
        IB(p).zf.reg_adjust.rotmap = zfreg_rot; IB(p).zf.local_adjust.rotmap = zfloc_rot; clear zf*_rot;
        DEM_diff.reg_adjust = IB(p).zo.reg_adjust.map - IB(p).zf.reg_adjust.rotmap;
        DEM_diff.local_adjust = IB(p).zo.local_adjust.map - IB(p).zf.local_adjust.rotmap;
        
        %filter the elevation difference maps to remove erroneous values
        DEM_diff.reg_adjust(DEM_diff.reg_adjust == 0) = NaN;
        DEM_diff.reg_adjust(isnan(IB(p).zo.reg_adjust.map)) = NaN; DEM_diff.reg_adjust(isnan(IB(p).zf.reg_adjust.rotmap)) = NaN;
        DEM_diff.local_adjust(DEM_diff.local_adjust == 0) = NaN;
        DEM_diff.local_adjust(isnan(IB(p).zo.local_adjust.map)) = NaN; DEM_diff.local_adjust(isnan(IB(p).zf.local_adjust.rotmap)) = NaN;
        disp('Extracted elevation differences from (regional & local) sea level-adjusted DEMs');
        
        %calculate the mean & median elevation within each polygon for comparison
        non_nans = ~isnan(IB(p).zo.local_adjust.map) & ~isnan(IB(p).zf.local_adjust.rotmap);
        IB(p).zo.reg_adjust.mean = nanmean(IB(p).zo.reg_adjust.map(non_nans));
        IB(p).zf.reg_adjust.mean = nanmean(IB(p).zf.reg_adjust.rotmap(non_nans));
        IB(p).zo.local_adjust.mean = nanmean(IB(p).zo.local_adjust.map(non_nans));
        IB(p).zf.local_adjust.mean = nanmean(IB(p).zf.local_adjust.rotmap(non_nans));
        
        %filter-out elevation change outliers & save the raw and filtered
        %sea level-adjusted data to the structure
        IB(p).dz.reg_adjust.map = DEM_diff.reg_adjust;
        IB(p).dz.local_adjust.map = DEM_diff.local_adjust;
        IB(p).dz.reg_adjust.mean = nanmean(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.mean = nanmean(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        IB(p).dz.reg_adjust.stdev = nanstd(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.stdev = nanstd(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        IB(p).dz.reg_adjust.median = nanmedian(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.median = nanmedian(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        IB(p).dz.reg_adjust.mad = mad(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)),1);
        IB(p).dz.local_adjust.mad = mad(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)),1);
        IB(p).dz.reg_adjust.max = max(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.max = max(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        IB(p).dz.reg_adjust.min = min(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.min = min(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        non_outliers = IB(p).dz.local_adjust.map<IB(p).dz.local_adjust.median+3*(1.4826*IB(p).dz.local_adjust.mad) & IB(p).dz.local_adjust.map>IB(p).dz.local_adjust.median-3*(1.4826*IB(p).dz.local_adjust.mad);
        IB(p).dz.reg_adjust.filtered_mean = nanmean(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_mean = nanmean(IB(p).dz.local_adjust.map(non_outliers));
        IB(p).dz.reg_adjust.filtered_stdev = nanstd(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_stdev = nanstd(IB(p).dz.local_adjust.map(non_outliers));
        IB(p).dz.reg_adjust.filtered_median = nanmedian(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_median = nanmedian(IB(p).dz.local_adjust.map(non_outliers));
        IB(p).dz.reg_adjust.filtered_mad = mad(IB(p).dz.reg_adjust.map(non_outliers & ~isnan(IB(p).dz.reg_adjust.map)),1);
        IB(p).dz.local_adjust.filtered_mad = mad(IB(p).dz.local_adjust.map(non_outliers & ~isnan(IB(p).dz.local_adjust.map)),1);
        IB(p).dz.reg_adjust.filtered_max = max(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_max = max(IB(p).dz.local_adjust.map(non_outliers));
        IB(p).dz.reg_adjust.filtered_min = min(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_min = min(IB(p).dz.local_adjust.map(non_outliers));
        
        %display select values
        disp(['Mean initial freeboard ',num2str(IB(p).zo.local_adjust.mean)]);
        disp(['Mean dz (locally-coregistered) ',num2str(IB(p).dz.local_adjust.filtered_mean),'m']);
        disp(['Mean dz (regionally-coregistered) ',num2str(IB(p).dz.reg_adjust.filtered_mean),'m']);
        
        %plot the elevation change map to determine whether the translation
        %and rotation estimates are reasonably accurate
        disp('Check the elevation change map: marginal areas of big change mean the mask is bad');
        figA=figure; imagesc(IB(p).xo,IB(p).yo,IB(p).dz.reg_adjust.map); axis xy equal; colormap jet; colorbar;
%         set(gca,'clim',[nanmedian(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)))-3*1.4826*mad(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)),1) nanmedian(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)))+3*1.4826*mad(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)),1)]);
        set(gca,'clim',[nanmean(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)))-2*nanstd(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map))) nanmean(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)))+2*nanstd(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)))]);
        prompt = 'Is the range of elevation change HUGE across the iceberg (i.e., tens of meters) (y/n)?'; str = input(prompt,'s');
        if strmatch(str,'n')==1
            close(figA);
            %advance the loop or recalculate sea level if melt rates are negative
            if p==1 && IB(p).dz.local_adjust.filtered_mean <=0 
                disp('Elevation change suggests ice accretion. Re-calculate elevation adjustment');
                clear IB;
                prompt = 'Should the regional, NOT local, fjord elevation offsets be used to coregister the DEMs (y/n)?';
                str = input(prompt,'s');
                if strmatch(str,'y')==1
                    sl_offset1 = nanmedian(min(fjordz_o)); sl_offset2 = nanmedian(min(fjordz_f));
                    A.z_local_adjust = A.z_reg_adjust; B.z_local_adjust = B.z_reg_adjust; 
                    flag = 1;
                else
                    disp('...then re-extract local sea level elevations from the DEMs');
                    [A,B,sl_offset1,sl_offset2,~,~,flag] = sea_level_adjust(DEM1.time,DEM2.time,IM1,IM2,dir_output,A,B,S,sl_early,sl_late,fjord_offset,fjordz_o,fjordz_f,iceberg_no,region_abbrev,cmax,DEMo_pos,DEMf_pos,imo_pos,imf_pos);
                end
                p = 1;
            else
                if p == 5
                    close all; break
                else
                    p = p + 1; close all;
                end
                drawnow;
                clear new_ang x_rot y_rot;
            end
        else
            disp('recalculate iceberg translation & rotation');
            if p ==1
                clear IB
            else
            ib = IB(1:p-1); clear IB; IB = ib; clear ib;
            end
            close(figA); p = p; 
        end
    elseif p > 1 && ~isempty(strmatch('Sl',w))
        %all iterations will have the same values b/c the translation &
        %rotation are pre-determined based on the polygon shape-files
        IB(p).xo = IB(1).xo; IB(p).yo = IB(1).yo;
        IB(p).zo.reg_adjust.map = IB(1).zo.reg_adjust.map;
        IB(p).zo.local_adjust.map = IB(1).zo.local_adjust.map;
        IB(p).xf = IB(1).xf; IB(p).yf = IB(1).yf;
        IB(p).zf.reg_adjust.map = IB(1).zf.reg_adjust.map;
        IB(p).zf.local_adjust.map = IB(1).zf.local_adjust.map;
        IB(p).zf.local_adjust.rotmap = IB(1).zf.local_adjust.rotmap;
%         IB(p).zo.local_adjust.rotmap = IB(1).zo.local_adjust.rotmap;
        IB(p).zo.local_adjust.mean = IB(1).zo.local_adjust.mean;
        IB(p).zf.local_adjust.mean = IB(1).zf.local_adjust.mean;
        
        %filter-out elevation change outliers & save the raw and filtered
        %sea level-adjusted data to the structure
        IB(p).dz.reg_adjust.map = IB(1).dz.reg_adjust.map;
        IB(p).dz.local_adjust.map = IB(1).dz.local_adjust.map;
        IB(p).dz.reg_adjust.mean = nanmean(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.mean = nanmean(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        IB(p).dz.reg_adjust.stdev = nanstd(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.stdev = nanstd(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        IB(p).dz.reg_adjust.median = nanmedian(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.median = nanmedian(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        IB(p).dz.reg_adjust.mad = mad(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)),1);
        IB(p).dz.local_adjust.mad = mad(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)),1);
        IB(p).dz.reg_adjust.max = max(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.max = max(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        IB(p).dz.reg_adjust.min = min(IB(p).dz.reg_adjust.map(~isnan(IB(p).dz.reg_adjust.map)));
        IB(p).dz.local_adjust.min = min(IB(p).dz.local_adjust.map(~isnan(IB(p).dz.local_adjust.map)));
        non_outliers = IB(p).dz.local_adjust.map<IB(p).dz.local_adjust.median+3*(1.4826*IB(p).dz.local_adjust.mad) & IB(p).dz.local_adjust.map>IB(p).dz.local_adjust.median-3*(1.4826*IB(p).dz.local_adjust.mad);
        IB(p).dz.reg_adjust.filtered_mean = nanmean(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_mean = nanmean(IB(p).dz.local_adjust.map(non_outliers));
        IB(p).dz.reg_adjust.filtered_stdev = nanstd(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_stdev = nanstd(IB(p).dz.local_adjust.map(non_outliers));
        IB(p).dz.reg_adjust.filtered_median = nanmedian(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_median = nanmedian(IB(p).dz.local_adjust.map(non_outliers));
        IB(p).dz.reg_adjust.filtered_mad = mad(IB(p).dz.reg_adjust.map(non_outliers & ~isnan(IB(p).dz.reg_adjust.map)),1);
        IB(p).dz.local_adjust.filtered_mad = mad(IB(p).dz.local_adjust.map(non_outliers & ~isnan(IB(p).dz.local_adjust.map)),1);
        IB(p).dz.reg_adjust.filtered_max = max(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_max = max(IB(p).dz.local_adjust.map(non_outliers));
        IB(p).dz.reg_adjust.filtered_min = min(IB(p).dz.reg_adjust.map(non_outliers));
        IB(p).dz.local_adjust.filtered_min = min(IB(p).dz.local_adjust.map(non_outliers));

        %decide if you need to re-do the rotation step or advance the loop
        if IB(p).dz.local_adjust.filtered_mean <0
            disp('');
            disp('Elevation change suggests ice accretion. Re-calculate iceberg rotation');
            p = p;
        else
            if p == 5; break; else p = p + 1; close all; end
        end
        clear new_ang x_rot y_rot median_dz;
    end
end
save([dir_output,'/',DEM1.time,'-',DEM2.time,'/','iceberg',iceberg_no,'_dz.mat'],'IB');
disp('iceberg info saved');

%plot the average elevation change map
for i = 1:length(IB)
    dzmap(:,:,i) = IB(i).dz.local_adjust.map;
    dzmean(i) = IB(i).dz.local_adjust.filtered_mean;
end
figure; set(gcf,'position',[50 50 800 500]);
subplot(1,2,1);
imagesc(IB(p).xo,IB(p).yo,IB(p).zo.local_adjust.map); axis xy equal; colormap jet; colorbar;
subplot(1,2,2);
imagesc(IB(p).xo,IB(p).yo,nanmean(dzmap,3)); axis xy equal; colormap jet; colorbar;
if nanmean(dzmean) < 0; flag = 1; end %use regional sea level correction if mean elevation change suggests accretion

%plot a histogram
disp('Done estimating elevation change. Plot a histogram of elevation change using data from all iterations.')
close all; drawnow;
reg_vals = []; local_vals = []; %brt_vals = [];
for i = 1:length(IB)
    non_NaNs = length(IB(i).dz.reg_adjust.map(~isnan(IB(i).dz.reg_adjust.map)));
    reg_vals(size(reg_vals,2)+1:size(reg_vals,2)+non_NaNs) = IB(i).dz.reg_adjust.map(~isnan(IB(i).dz.reg_adjust.map));
    clear non_NaNs;
    non_NaNs = length(IB(i).dz.local_adjust.map(~isnan(IB(i).dz.local_adjust.map)));
    local_vals(size(local_vals,2)+1:size(local_vals,2)+non_NaNs) = IB(i).dz.local_adjust.map(~isnan(IB(i).dz.local_adjust.map));
    clear non_NaNs;
end
mean_regvals = nanmean(reg_vals); std_regvals = nanstd(reg_vals);
median_regvals = nanmedian(reg_vals); mad_regvals = mad(reg_vals,1);
reg_outliers = reg_vals>median_regvals+3*1.4826*mad_regvals | reg_vals<median_regvals-3*1.4826*mad_regvals;
reg_vals(reg_outliers) = NaN; min_regvals = min(reg_vals); max_regvals = max(reg_vals);
mean_locvals = nanmean(local_vals); std_locvals = nanstd(local_vals);
median_locvals = nanmedian(local_vals); mad_locvals = mad(local_vals,1);
local_outliers = local_vals>median_locvals+3*1.4826*mad_locvals | local_vals<median_locvals-3*1.4826*mad_locvals;
local_vals(local_outliers) = NaN; min_locvals = min(local_vals); max_locvals = max(local_vals);
figure; set(gcf,'position',[50 50 400 400]);
hist(local_vals,[min_locvals:(max_locvals-min_locvals)/50:max_locvals]);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k','EdgeColor','k');
hold on;
pl(1) = plot([mean_locvals mean_locvals],[0 max(get(gca,'ylim'))],'-r','linewidth',2);
pl(2) = plot([median_locvals median_locvals],[0 max(get(gca,'ylim'))],'-b','linewidth',2);
pl(3) = plot([nanmean(local_vals) nanmean(local_vals)],[0 max(get(gca,'ylim'))],'--r','linewidth',2);
pl(4) = plot([nanmedian(local_vals) nanmedian(local_vals)],[0 max(get(gca,'ylim'))],'--b','linewidth',2);
legend(pl,'unfiltered mean','unfiltered median','filtered mean','filtered median');
title('Local sea level-adjusted \Deltaz');
figure; set(gcf,'position',[450 50 400 400]);
hist(reg_vals,[min_regvals:(max_regvals-min_regvals)/50:max_regvals]);
h = findobj(gca,'Type','patch');
set(h,'FaceColor','k','EdgeColor','k');
hold on;
pl(1) = plot([mean_regvals mean_regvals],[0 max(get(gca,'ylim'))],'-r','linewidth',2);
pl(2) = plot([median_regvals median_regvals],[0 max(get(gca,'ylim'))],'-b','linewidth',2);
pl(3) = plot([nanmean(reg_vals) nanmean(reg_vals)],[0 max(get(gca,'ylim'))],'--r','linewidth',2);
pl(4) = plot([nanmedian(reg_vals) nanmedian(reg_vals)],[0 max(get(gca,'ylim'))],'--b','linewidth',2);
legend(pl,'unfiltered mean','unfiltered median','filtered mean','filtered median');
title('Regional sea level-adjusted \Deltaz');

if (sl_early.quality + sl_late.quality) == 2 || flag == 0 %use local sea level-adjusted elevations
    %display the mean elevation change & st. dev. used in volume flux calculations
    disp('Elevation change stats from local sea level-adjusted DEMs:');
    disp(['Mean filtered elevation change: ',num2str(nanmean(local_vals)),'m']);
    disp(['Variability about the mean (std): ',num2str(nanstd(local_vals)),'m']);
    
    %save the elevation change data
    dz.best_source='local_sealevel';
    dz.local.mean = nanmean(local_vals); dz.local.stdev = nanstd(local_vals);
    dz.local.median = nanmedian(local_vals); dz.local.mad = mad(local_vals,1);
    dz.regional.mean = nanmean(reg_vals); dz.regional.stdev = nanstd(reg_vals);
    dz.regional.median = nanmedian(reg_vals); dz.regional.mad = mad(reg_vals,1);
    save([dir_output,'/',DEM1.time,'-',DEM2.time,'/','iceberg',num2str(iceberg_no),'_dz.mat'],'IB','dz');
    disp('iceberg info with elevation change (dz) saved');    

else %use regional sea level-adjusted elevations
    %display the mean elevation change & st. dev. used in volume flux calculations
    disp('Elevation change stats from regional sea level-adjusted DEMs:');
    disp(['Mean filtered elevation change: ',num2str(nanmean(reg_vals)),'m']);
    disp(['Variability about the mean (std): ',num2str(nanstd(reg_vals)),'m']);
    
    %save the elevation change data
    dz.best_source='regional_sealevel';
    dz.local.mean = nanmean(local_vals); dz.local.stdev = nanstd(local_vals);
    dz.local.median = nanmedian(local_vals); dz.local.mad = mad(local_vals,1);
    dz.regional.mean = nanmean(reg_vals); dz.regional.stdev = nanstd(reg_vals);
    dz.regional.median = nanmedian(reg_vals); dz.regional.mad = mad(reg_vals,1);
    save([dir_output,'/',DEM1.time,'-',DEM2.time,'/','iceberg',num2str(iceberg_no),'_dz.mat'],'IB','dz');
    disp('iceberg info with elevation change (dz) saved');
end
    
%export the average mask for the iceberg in the later DEM
for i = 1:length(IB)
    xf(i,:) = IB(i).vertices.xf; yf(i,:) = IB(i).vertices.yf;
end
mean_xf = nanmean(xf); mean_yf = nanmean(yf);
cd([dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_rois/']);
clear S
S.Geometry = 'Polygon';
S.BoundingBox = [min(mean_xf) min(mean_yf); max(mean_xf) max(mean_yf)];
S.X = double(mean_xf);
S.Y = double(mean_yf);
S.Name = ['iceberg',num2str(iceberg_no)];
% shapefile_name = [region_abbrev,'_',num2str(DEM2.time),'_iceberg',num2str(iceberg_no)];
shapefile_name = [region_abbrev,'_',num2str(DEM2.time),'_iceberg',num2str(iceberg_no)];
shapewrite(S,shapefile_name);
copyfile([dir_code,PSprojfile],[dir_output,'/',DEM1.time,'-',DEM2.time,'/iceberg_rois/',shapefile_name,'.prj']);

%advance to the next iceberg
disp('Click anywhere on the last histogram figure to advance to the next iceberg...');
wait = waitforbuttonpress;
if wait == 0
    disp(['Saved iceberg #',num2str(iceberg_no),' elevation change data']);
end


end
 
