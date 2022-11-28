function [DEM1,DEM2] = calculate_DEM_offset(varargin)
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
%           dir_TMD         directory path to the TMD folder ('*/TMD/')
%           dir_output      directory where all output files will be placed
%           fjord_lat       central latitude coordinate of the fjord
%           fjord_lon       central longitude coordinate of the fjord
%
% OUTPUTS:  DEM1            DEM1 structure variable with new fields
%           DEM2            DEM2 structure variable with new fields
%
% Calls the following external functions if prompted:
%   - tmd_tide_pred.m

%assign variable names & specify options based on inputs
DEM1 = varargin{1}; DEM2 = varargin{2};
IM1 = varargin{3}; IM2 = varargin{4}; 
dir_output = varargin{5};
fjord_lon = varargin{6}; fjord_lat = varargin{7};
if nargin == 8
    dir_TMD = varargin{8};
end

%identify over-lapping bedrock regions
x_overlap = find(DEM2.x > min(DEM1.x) & DEM2.x < max(DEM1.x));
y_overlap = find(DEM2.y > min(DEM1.y) & DEM2.y < max(DEM1.y));
figure1 = figure; elev_cmap = colormap(jet(10001)); elev_cmap(1,:) = [1 1 1];
DEM1.z_masked = DEM1.z_elpsd_adjust; DEM1.z_masked(DEM1.mask==0)=NaN;
imagesc(DEM1.x,DEM1.y,DEM1.z_masked); hold on;
set(gca,'ydir','normal'); colormap(gca,elev_cmap); set(gca,'clim',[min(min(DEM1.z_masked)) min(min(DEM1.z_masked))+300]);
set(figure1,'position',[50 50 800 800]);
disp('To draw a polygon, wait for crosshairs to appear on the figure. Be patient!');
disp('    Click on the polygon vertices one at a time');
disp('    Adjust your vertices after completing the polygon by clicking and dragging.');
disp('    When you are satisfied with your polygon, right-click on the figure and select "Create Mask"');
disp('Draw a polygon around the map extent containing elevations in the earlier DEM');
Z_extent = roipoly;
clf(figure1);
disp('Draw a polygon around the map extent containing elevations in the later DEM');
DEM2.z_masked = DEM2.z_elpsd_adjust; DEM2.z_masked(DEM2.mask==0)=NaN;
imagesc(DEM2.x,DEM2.y,DEM2.z_masked); hold on;
set(gca,'ydir','normal'); colormap(gca,elev_cmap); set(gca,'clim',[min(min(DEM2.z_masked)) min(min(DEM2.z_masked))+300]);
Y_extent = roipoly;
[Yx,Yy] = meshgrid(DEM2.x,DEM2.y); [Zx,Zy] = meshgrid(DEM1.x,DEM1.y);
Z_extent_interp = interp2(Zx,Zy,double(Z_extent),Yx,Yy);
overlap = Z_extent_interp + double(Y_extent);
overlap(overlap<2) = 0; overlap(overlap==2) = 1;
z_overlap = overlap.*DEM2.z_masked;
close(figure1);

%plot both the images for help identifying bedrock regions & thin sea ice in the fjord
I1.x = IM1.x(1:2:end); I1.y = IM1.y(1:2:end);
I1.z = IM1.z(1:2:end,1:2:end);
I1.z_masked = IM1.z(1:2:end,1:2:end);
I2.x = IM2.x(1:2:end); I2.y = IM2.y(1:2:end);
I2.z = IM2.z(1:2:end,1:2:end);
I2.z_masked = IM2.z(1:2:end,1:2:end);
figureA = figure; set(figureA,'position',[850 50 400 900]);
sub1 = subplot(2,1,1);
imagesc(I1.x,I1.y,I1.z_masked); axis xy equal; colormap gray; hold on;
sub2 = subplot(2,1,2);
imagesc(I2.x,I2.y,I2.z_masked); axis xy equal; colormap gray; hold on;

%identify fjord and bedrock regions of overlap
if isnan(nanmean(z_overlap(z_overlap>0))) || isempty(y_overlap) || isempty(x_overlap) || double(size(y_overlap,2)./size(DEM2.z_masked,1))<0.1 || double(size(x_overlap,2)./size(DEM2.z_masked,2))<0.1;
    disp('No bedrock overlap: use sea level to adjust fjord elevations');
    
    % adjust sea level elevations using ice-free fjord elevations
    %Early image
    disp('Draw a polygon around the ice-free portion of the fjord in each DEM');
    figure1 = figure;
    imagesc(DEM1.x,DEM1.y,DEM1.z_masked); hold on;
    colormap(gca,elev_cmap);
    set(gca,'ydir','normal'); set(gca,'clim',[min(min(DEM1.z_masked)) min(min(DEM1.z_masked))+300]);
    set(figure1,'position',[50 50 800 800]);
%     figure1b = figure; set(gcf,'position',[850 50 800 800]);
%     imagesc(IM.x,IM.y,IM.z); colormap gray;
%     set(gca,'ydir','normal');
    disp('Early DEM');
    figure(figure1);
    [fjord_mask_early,xv,yv] = roipoly;
    fjord.dem1.z = fjord_mask_early.*DEM1.z_masked; fjord.dem1.z(fjord.dem1.z == 0) = NaN;
    fjord.dem1.x = nearestneighbour(xv',DEM1.x); fjord.dem1.y = nearestneighbour(yv',DEM1.y);
    clear xv yv; clear IM; 
    %Late image
    disp('Later DEM');
    figure2 = figure;
    imagesc(DEM2.x,DEM2.y,DEM2.z_masked); hold on;
    colormap(gca,elev_cmap);
    set(gca,'ydir','normal'); set(gca,'clim',[min(min(DEM2.z_masked)) min(min(DEM2.z_masked))+300]);
    set(figure2,'position',[50 50 800 800]);
%     figure2b = figure; set(gcf,'position',[850 50 800 800]);
%     imagesc(IM.x,IM.y,IM.z); colormap gray;
%     set(gca,'ydir','normal');
    figure(figure2);
    [fjord_mask_late,xv,yv] = roipoly;
    fjord.dem2.z = fjord_mask_late.*DEM2.z_masked; fjord.dem2.z(fjord.dem2.z == 0) = NaN;
    fjord.dem2.x = nearestneighbour(xv',DEM2.x); fjord.dem2.y = nearestneighbour(yv',DEM2.y);
    clear xv yv; clear IM;
    close(figure1); close(figure2);
    
    %rename variables
    z1 = fjord.dem1.z(~isnan(fjord.dem1.z)); z2 = fjord.dem2.z(~isnan(fjord.dem2.z));
    
    %initial stats
    z1_sl = min(z1); z2_sl = min(z2);
    z1_med = median(z1); z1_mad = mad(z1);
    z2_med = median(z2); z2_mad = mad(z2);
    
    %filter
    z1(z1>z1_med+1.48*z1_mad | z1<z1_med-1.48*z1_mad) = NaN;
    z2(z2>z2_med+1.48*z2_mad | z2<z2_med-1.48*z2_mad) = NaN;
    
    %recalculate stats
    z1_med = nanmedian(z1); z1_mad = mad(z1);
    z2_med = nanmedian(z2); z2_mad = mad(z2);
    z1_mean = nanmean(z1);
    z2_mean = nanmean(z2);
    fjord.dem_diff.mean = z1_mean-z2_mean;
    fjord.dem_diff.median = z1_med-z2_med;
    fjord.dem_diff.sl_adjust = fjord.dem_diff.median;
    disp(['Sea level offset (sea level adjusted): ',num2str(fjord.dem_diff.median)]);
    
    % calculate the tidal change between acquistion dates if desired
    if nargin == 8
        disp('Extract tidal change between acquisition dates');
        early_time = datenum(DEM1.time,'yyyymmddHHMMSS');
        late_time = datenum(DEM2.time,'yyyymmddHHMMSS');
        SD_time = [early_time:(late_time-early_time)/10:late_time];
        cd([dir_TMD,'DATA/']);
        [tidal_z,~] = tmd_tide_pred('Model_CATS2008',SD_time,fjord_lat,fjord_lon,'z'); %change in tide at the fjord mouth
        DEM1.date = SD_time(1); DEM2.date = SD_time(end);
        fjord.tidal_change = tidal_z(1)-tidal_z(end); %tidal change between image acquisition times
    end
    
    %save the fjord offset data
    save([dir_output,'fjord_offset_',DEM1.time,'-',DEM2.time,'.mat'],'fjord');

    %add the image date to each structure
    early_time = datenum(num2str(DEM1.time),'yyyymmddHHMMSS');
    late_time = datenum(num2str(DEM2.time),'yyyymmddHHMMSS');
    DEM1.date = early_time; DEM2.date = late_time;
    
    %adjust sea level (coregister DEMs)
    DEM2.z_masked_sl_adjust = DEM2.z_masked + fjord.dem_diff.sl_adjust;
    close all;
    
    %save the DEMs
    disp('Resaving DEMs with the sea level adjustment info added to their structures');
    DEM=DEM1; save([dir_output,DEM1.filename,'.mat'],'DEM','-append'); clear DEM; 
    DEM=DEM2; save([dir_output,DEM2.filename,'.mat'],'DEM','-append'); clear DEM;
    
else
    %find the region of DEM overlap
    minZx = min(DEM1.x); maxZx = max(DEM1.x); minZy = min(DEM1.y); maxZy = max(DEM1.y);
    minYx = min(DEM2.x); maxYx = max(DEM2.x); minYy = min(DEM2.y); maxYy = max(DEM2.y);
    minx = [minZx minYx]; maxx = [maxZx maxYx];
    miny = [minZy minYy]; maxy = [maxZy maxYy];
    x_bounds = [max(minx) min(maxx)]; y_bounds = [max(miny) min(maxy)];
    
    % adjust sea level elevations using bedrock offset & tide change and ice-free fjord elevations
    figure1 = figure; set(gcf,'position',[50 50 800 800]);
    overlap_elevs = overlap.*DEM2.z_masked; overlap_elevs(isnan(DEM2.z_masked)) = NaN; overlap_elevs(overlap==0) = NaN;
    imagesc(DEM2.x,DEM2.y,overlap_elevs);
    set(gca,'ydir','normal'); colormap(gca,elev_cmap); 
    set(gca,'clim',[min(min(DEM2.z_masked)) min(min(DEM2.z_masked))+300]);
    set(gca,'xlim',x_bounds,'ylim',y_bounds);
%     figure1b = figure; set(gcf,'position',[850 50 800 800]);
%     imagesc(IM.x,IM.y,IM.z); colormap gray;
%     set(gca,'ydir','normal'); 
    figure(figureA); 
    subplot(sub1); set(gca,'xlim',x_bounds,'ylim',y_bounds);
    subplot(sub2); set(gca,'xlim',x_bounds,'ylim',y_bounds);
    
    
    %confirm that there is bedrock in the overlapping portion of the images
    prompt = 'Is there bedrock in the overlapping DEM (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        
        %Region 1 mask
        disp('Draw 2 polygons over bedrock, preferably on opposite sides of the fjord');
        disp('Draw a polygon around the 1st overlapping bedrock region');
        figure(figure1);
        [bedmask_reg1_Y,xv,yv] = roipoly;
        bedrock_reg1_Y = bedmask_reg1_Y.*DEM2.z_masked; bedrock_reg1_Y(bedrock_reg1_Y == 0) = NaN;
        x1_Y = nearestneighbour(xv',DEM2.x); y1_Y = nearestneighbour(yv',DEM2.y);
        x1_Z = nearestneighbour(xv',DEM1.x); y1_Z = nearestneighbour(yv',DEM1.y);
        bedmask_reg1_Z = roipoly(DEM1.z_masked,x1_Z,y1_Z);
        bedrock_reg1_Z = bedmask_reg1_Z.*DEM1.z_masked; bedrock_reg1_Z(bedrock_reg1_Z == 0) = NaN;
        clear xv yv;
        %Region 2 mask
        disp('Draw a polygon around the 2nd overlapping bedrock region');
        [bedmask_reg2_Y,xv,yv] = roipoly;
        bedrock_reg2_Y = bedmask_reg2_Y.*DEM2.z_masked; bedrock_reg2_Y(bedrock_reg2_Y == 0) = NaN;
        x2_Y = nearestneighbour(xv',DEM2.x); y2_Y = nearestneighbour(yv',DEM2.y);
        x2_Z = nearestneighbour(xv',DEM1.x); y2_Z = nearestneighbour(yv',DEM1.y);
        bedmask_reg2_Z = roipoly(DEM1.z_masked,x2_Z,y2_Z);
        bedrock_reg2_Z = bedmask_reg2_Z.*DEM1.z_masked; bedrock_reg2_Z(bedrock_reg2_Z == 0) = NaN;
        
        %Region 1
        bedrock.reg1.dem1.x = DEM1.x(1,min(x1_Z):max(x1_Z)); bedrock.reg1.dem1.y = DEM1.y(1,min(y1_Z):max(y1_Z));
        bedrock.reg1.dem1.z = bedrock_reg1_Z(min(y1_Z):max(y1_Z),min(x1_Z):max(x1_Z));
        bedrock.reg1.dem2.x = DEM2.x(1,min(x1_Y):max(x1_Y)); bedrock.reg1.dem2.y = DEM2.y(1,min(y1_Y):max(y1_Y));
        bedrock.reg1.dem2.z = bedrock_reg1_Y(min(y1_Y):max(y1_Y),min(x1_Y):max(x1_Y));
        %Region 2
        bedrock.reg2.dem1.x = DEM1.x(1,min(x2_Z):max(x2_Z)); bedrock.reg2.dem1.y = DEM1.y(1,min(y2_Z):max(y2_Z));
        bedrock.reg2.dem1.z = bedrock_reg2_Z(min(y2_Z):max(y2_Z),min(x2_Z):max(x2_Z));
        bedrock.reg2.dem2.x = DEM2.x(1,min(x2_Y):max(x2_Y)); bedrock.reg2.dem2.y = DEM2.y(1,min(y2_Y):max(y2_Y));
        bedrock.reg2.dem2.z = bedrock_reg2_Y(min(y2_Y):max(y2_Y),min(x2_Y):max(x2_Y));
        
        
        %interpolate to the same grid
        %Region 1
        [x_grid_dem1,y_grid_dem1] = meshgrid(bedrock.reg1.dem1.x,bedrock.reg1.dem1.y);
        [x_grid_dem2,y_grid_dem2] = meshgrid(bedrock.reg1.dem2.x,bedrock.reg1.dem2.y);
        z = interp2(x_grid_dem2,y_grid_dem2,bedrock.reg1.dem2.z,x_grid_dem1,y_grid_dem1);
        bedrock.reg1.dem2.z_interp = z;
        clear x_grid_dem1 y_grid_dem1 x_grid_dem2 y_grid_dem2 z;
        %Region 2
        [x_grid_dem1,y_grid_dem1] = meshgrid(bedrock.reg2.dem1.x,bedrock.reg2.dem1.y);
        [x_grid_dem2,y_grid_dem2] = meshgrid(bedrock.reg2.dem2.x,bedrock.reg2.dem2.y);
        z = interp2(x_grid_dem2,y_grid_dem2,bedrock.reg2.dem2.z,x_grid_dem1,y_grid_dem1);
        bedrock.reg2.dem2.z_interp = z;
        clear x_grid_dem1 y_grid_dem1 x_grid_dem2 y_grid_dem2 z;
        
        %crop values where data intersect an image boundary
        %Region 1
        border = bedrock.reg1.dem1.z == 0 | bedrock.reg1.dem2.z_interp == 0;
        bedrock.reg1.dem2.z_interp(border) = 0; bedrock.reg1.dem1.z(border) = 0;
        nans = isnan(bedrock.reg1.dem1.z) | isnan(bedrock.reg1.dem2.z_interp);
        bedrock.reg1.dem2.z_interp(nans) = 0; bedrock.reg1.dem1.z(nans) = 0;
        clear border nans;
        %Region 2
        border = bedrock.reg2.dem1.z == 0 | bedrock.reg2.dem2.z_interp == 0;
        bedrock.reg2.dem2.z_interp(border) = 0; bedrock.reg2.dem1.z(border) = 0;
        nans = isnan(bedrock.reg2.dem1.z) | isnan(bedrock.reg2.dem2.z_interp);
        bedrock.reg2.dem2.z_interp(nans) = 0; bedrock.reg2.dem1.z(nans) = 0;
        clear border nans;
        
        %plot the subsetted regions
        %Region 1
        clf(figure1); %clf(figure1b);
        subset1_elevmax = max(max(bedrock.reg1.dem1.z)); subset2_elevmax = max(max(bedrock.reg2.dem1.z));
        figure(figure1);
        imagesc(bedrock.reg1.dem1.x,bedrock.reg1.dem1.y,bedrock.reg1.dem1.z)
        set(gca,'ydir','normal','clim',[0 min([subset1_elevmax 1000])]); colorbar;
        set(figure1,'position',[50 550 500 500]);
        figure2 = figure;
        imagesc(bedrock.reg1.dem1.x,bedrock.reg1.dem1.y,bedrock.reg1.dem2.z_interp)
        set(gca,'ydir','normal','clim',[0 min([subset1_elevmax 1000])]); colorbar;
        set(figure2,'position',[550 550 500 500]);
        %Region 2
        figure3 = figure;
        imagesc(bedrock.reg2.dem1.x,bedrock.reg2.dem1.y,bedrock.reg2.dem1.z)
        set(gca,'ydir','normal','clim',[0 min([subset2_elevmax 1000])]); colorbar;
        set(figure3,'position',[50 0 500 500]);
        figure4 = figure;
        imagesc(bedrock.reg2.dem1.x,bedrock.reg2.dem1.y,bedrock.reg2.dem2.z_interp)
        set(gca,'ydir','normal','clim',[0 min([subset2_elevmax 1000])]); colorbar;
        set(figure4,'position',[550 0 500 500]);
        
        %difference the subsets
        %Region 1
        bedrock.reg1.dem_diff.map = bedrock.reg1.dem1.z-bedrock.reg1.dem2.z_interp;
        bedrock.reg1.dem_diff.map(bedrock.reg1.dem1.z == 0 & bedrock.reg1.dem2.z_interp == 0) = NaN;
        %Region 2
        bedrock.reg2.dem_diff.map = bedrock.reg2.dem1.z-bedrock.reg2.dem2.z_interp;
        bedrock.reg2.dem_diff.map(bedrock.reg2.dem1.z == 0 & bedrock.reg2.dem2.z_interp == 0) = NaN;
        
        %plot the DEM difference maps to look for spatial trends in bedrock offset
        %Region 1
        reg1_mean = nanmean(bedrock.reg1.dem_diff.map(~isnan(bedrock.reg1.dem_diff.map)));
        reg2_mean = nanmean(bedrock.reg2.dem_diff.map(~isnan(bedrock.reg2.dem_diff.map)));
        avg_offset = (reg1_mean + reg2_mean)/2;
        figure5 = figure;
        imagesc(bedrock.reg1.dem1.x,bedrock.reg1.dem1.y,bedrock.reg1.dem_diff.map);
        set(gca,'ydir','normal'); colorbar;
        set(gca,'clim',[avg_offset-5 avg_offset+5]);
        set(figure5,'position',[1050 550 500 500]);
        %Region 2
        figure6 = figure;
        imagesc(bedrock.reg2.dem1.x,bedrock.reg2.dem1.y,bedrock.reg2.dem_diff.map);
        set(gca,'ydir','normal'); colorbar;
        set(gca,'clim',[avg_offset-5 avg_offset+5]);
        set(figure6,'position',[1050 0 500 500]);
        
        %calculate statistics for the bedrock offset
        %Region 1
        bedrock.reg1.dem_diff.mean = nanmean(bedrock.reg1.dem_diff.map(~isnan(bedrock.reg1.dem_diff.map)));
        bedrock.reg1.dem_diff.median = nanmedian(bedrock.reg1.dem_diff.map(~isnan(bedrock.reg1.dem_diff.map)));
        bedrock.reg1.dem_diff.stdev = (nansum(nanstd(bedrock.reg1.dem_diff.map).^2)/size(bedrock.reg1.dem_diff.map,2)).^(1/2);
        % disp(num2str(bedrock.reg1.dem_diff.mean));
        %Region 2
        bedrock.reg2.dem_diff.mean = nanmean(bedrock.reg2.dem_diff.map(~isnan(bedrock.reg2.dem_diff.map)));
        bedrock.reg2.dem_diff.median = nanmedian(bedrock.reg2.dem_diff.map(~isnan(bedrock.reg2.dem_diff.map)));
        bedrock.reg2.dem_diff.stdev = (nansum(nanstd(bedrock.reg2.dem_diff.map).^2)/size(bedrock.reg2.dem_diff.map,2)).^(1/2);
        % disp(num2str(bedrock.reg2.dem_diff.mean));
        
        %filter-out outliers (mean +/- 2*stdev) & re-calculate offset stats
        %Region 1
        reg1_filter = bedrock.reg1.dem_diff.map >= bedrock.reg1.dem_diff.mean+2*bedrock.reg1.dem_diff.stdev | bedrock.reg1.dem_diff.map <= bedrock.reg1.dem_diff.mean-2*bedrock.reg1.dem_diff.stdev;
        bedrock.reg1.dem_diff.map(reg1_filter) = NaN;
        bedrock.reg1.dem_diff.mean = nanmean(bedrock.reg1.dem_diff.map(~isnan(bedrock.reg1.dem_diff.map)));
        bedrock.reg1.dem_diff.median = nanmedian(bedrock.reg1.dem_diff.map(~isnan(bedrock.reg1.dem_diff.map)));
        bedrock.reg1.dem_diff.stdev = (nansum(nanstd(bedrock.reg1.dem_diff.map).^2)/size(bedrock.reg1.dem_diff.map,2)).^(1/2);
        disp(['Region 1 bedrock offset: ',num2str(bedrock.reg1.dem_diff.mean), 'm']);
        %Region 2
        reg2_filter = bedrock.reg2.dem_diff.map >= bedrock.reg2.dem_diff.mean+2*bedrock.reg2.dem_diff.stdev | bedrock.reg2.dem_diff.map <= bedrock.reg2.dem_diff.mean-2*bedrock.reg2.dem_diff.stdev;
        bedrock.reg2.dem_diff.map(reg2_filter) = NaN;
        bedrock.reg2.dem_diff.mean = nanmean(bedrock.reg2.dem_diff.map(~isnan(bedrock.reg2.dem_diff.map)));
        bedrock.reg2.dem_diff.median = nanmedian(bedrock.reg2.dem_diff.map(~isnan(bedrock.reg2.dem_diff.map)));
        bedrock.reg2.dem_diff.stdev = (nansum(nanstd(bedrock.reg2.dem_diff.map).^2)/size(bedrock.reg2.dem_diff.map,2)).^(1/2);
        disp(['Region 2 bedrock offset: ',num2str(bedrock.reg2.dem_diff.mean), 'm']);
        
        %plot filtered data
        avg_offset = (bedrock.reg1.dem_diff.mean + bedrock.reg2.dem_diff.mean)/2;
        figure(figure5);
        imagesc(bedrock.reg1.dem1.x,bedrock.reg1.dem1.y,bedrock.reg1.dem_diff.map);
        set(gca,'ydir','normal'); colorbar;
        set(gca,'clim',[avg_offset-5 avg_offset+5]);
        %Region 2
        figure(figure6);
        h = imagesc(bedrock.reg2.dem1.x,bedrock.reg2.dem1.y,bedrock.reg2.dem_diff.map);
        set(gca,'ydir','normal'); colorbar;
        set(gca,'clim',[avg_offset-5 avg_offset+5]);
        
        %save the bedrock comparison data
        save([dir_output,'bedrock_offset_',DEM1.time,'-',DEM2.time,'.mat'],'bedrock');
        
        % compare ice-free elevations in the fjord
        disp('Check figures then close the bottom right figure to continue');
        waitfor(h); 
        close(figure1); close(figure2); close(figure3); close(figure4); close(figure5); drawnow;
        figure1 = figure;
        set(figure1,'position',[50 50 800 800]);
        imagesc(DEM2.x,DEM2.y,overlap.*DEM2.z_masked); hold on;
        set(gca,'ydir','normal'); set(gca,'clim',[0,500]);
        set(gca,'xlim',x_bounds,'ylim',y_bounds);
%         figure1b = figure; set(gcf,'position',[850 50 800 800]);
%         imagesc(IM.x,IM.y,IM.z); colormap gray;
%         set(gca,'ydir','normal'); set(gca,'xlim',x_bounds,'ylim',y_bounds);
        clear xv yv;
        disp('Draw polygon around the overlapping portion of the fjord with open water');
        figure(figure1);
        [fjord_mask,xv,yv] = roipoly;
        fjord_mask = double(fjord_mask);
        x1 = nearestneighbour(xv',DEM1.x); y1 = nearestneighbour(yv',DEM1.y);
        fjord.dem1.x = DEM1.x(min(x1):max(x1)); fjord.dem1.y = DEM1.y(min(y1):max(y1));
        x2 = nearestneighbour(xv',DEM2.x); y2 = nearestneighbour(yv',DEM2.y);
        fjord.dem2.x = DEM2.x(min(x2):max(x2)); fjord.dem2.y = DEM2.y(min(y2):max(y2));
        z2 = fjord_mask.*DEM2.z_masked; z2(z2 == 0) = NaN;
        fjord.dem2.z = z2(min(y2):max(y2),min(x2):max(x2));
        fjord.dem1.z = DEM1.z_masked(min(y1):max(y1),min(x1):max(x1));
        fjord_mask_crop = fjord_mask(min(y2):max(y2),min(x2):max(x2));
        clear xv yv;
        
        %interpolate to the same grid
        [x_grid_dem1,y_grid_dem1] = meshgrid(fjord.dem1.x,fjord.dem1.y);
        [x_grid_dem2,y_grid_dem2] = meshgrid(fjord.dem2.x,fjord.dem2.y);
        fjord.dem1.z_interp = interp2(x_grid_dem1,y_grid_dem1,fjord.dem1.z,x_grid_dem2,y_grid_dem2);
        fjord.dem1.z_interp = fjord_mask_crop.*fjord.dem1.z_interp; fjord.dem1.z_interp(fjord.dem1.z_interp == 0) = NaN;
        
        %plot the fjord data
        close all; drawnow;
        cmin = nanmedian(min(fjord.dem2.z(fjord.dem2.z>0)));
        figure7 = figure;
        imagesc(fjord.dem2.x,fjord.dem2.y,fjord.dem1.z_interp);
        set(gca,'ydir','normal'); colorbar;
        set(gca,'clim',[cmin, cmin+60]);
        set(figure7,'position',[50 550 500 500]);
        figure8 = figure;
        imagesc(fjord.dem2.x,fjord.dem2.y,fjord.dem2.z);
        set(gca,'ydir','normal'); colorbar;
        set(gca,'clim',[cmin, cmin+60]);
        set(figure8,'position',[550 550 500 500]);
        
        %calculate statistics for the fjord offset
        fjord.dem_diff.map = fjord.dem1.z_interp - fjord.dem2.z;
        fjord.dem_diff.mean = nanmean(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        fjord.dem_diff.stdev = nanstd(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        fjord.dem_diff.median = nanmedian(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        fjord.dem_diff.mad = mad(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        
        %filter the fjord min elevation data
        fjord_filter = fjord.dem_diff.map >= fjord.dem_diff.median+1.48*fjord.dem_diff.mad | fjord.dem_diff.map <= fjord.dem_diff.median-1.48*fjord.dem_diff.mad;
        fjord.dem_diff.map(fjord_filter) = NaN;
        fjord.dem_diff.mean = nanmean(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        fjord.dem_diff.stdev = nanstd(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        fjord.dem_diff.median = nanmedian(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        fjord.dem_diff.mad = mad(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        disp(['Sea level offset: ',num2str(fjord.dem_diff.mean)]);
        
        %save the fjord offset data
        save([dir_output,'fjord_offset_',DEM1.time,'-',DEM2.time,'.mat'],'fjord');
        
        % calculate the sea level adjustment using bedrock to adjust DEM
        % bias & the tidal model is desired
        if nargin == 8
            disp('Extract tidal change between acquisition dates');
            early_time = datenum(DEM1.time,'yyyymmddHHMMSS');
            late_time = datenum(DEM2.time,'yyyymmddHHMMSS');
            SD_time = [early_time:(late_time-early_time)/10:late_time];
            cd([dir_TMD,'DATA/']);
            [tidal_z,~] = tmd_tide_pred('Model_CATS2008',SD_time,fjord_lat,fjord_lon,'z'); %change in tide at the fjord mouth
            DEM1.date = SD_time(1); DEM2.date = SD_time(end);
            fjord.tidal_change = tidal_z(1)-tidal_z(end); %tidal change between image acquisition times
            if isnan(fjord.tidal_change)
                fjord.tidal_change = 0;
            end
            
            %coregister
            mean_bedrock_offset = (bedrock.reg1.dem_diff.mean+bedrock.reg2.dem_diff.mean)/2;
            fjord.dem_diff.br_tide_adjust = mean_bedrock_offset + fjord.tidal_change;
            DEM2.z_masked_br_tide_adjust = DEM2.z_masked + fjord.dem_diff.br_tide_adjust;
            disp(['Sea level offset (bedrock & tide adjusted): ',num2str(fjord.dem_diff.br_tide_adjust), 'm']);
        end
        
        %coregister using the fjord minimum elevations
        close all;
        fjord.dem_diff.sl_adjust = fjord.dem_diff.median;
        DEM2.z_masked_sl_adjust = DEM2.z_masked + fjord.dem_diff.sl_adjust;
        disp(['Sea level offset (sea level adjusted): ',num2str(fjord.dem_diff.sl_adjust), 'm']);
        save([dir_output,'fjord_offset_',DEM1.time,'-',DEM2.time,'.mat'],'fjord');
        
    else %no bedrock in overlapping portion of DEMs
        
        % compare ice-free elevations in the fjord
        close(figure1);
        figure1 = figure;
        set(figure1,'position',[50 50 800 800]);
        overlap_elevs = overlap.*DEM2.z_masked; overlap_elevs(isnan(DEM2.z_masked)) = NaN; overlap_elevs(overlap==0) = NaN;
        imagesc(DEM2.x,DEM2.y,overlap_elevs);
        set(gca,'ydir','normal'); colormap(gca,elev_cmap); 
        set(gca,'clim',[min(min(DEM2.z_masked)) min(min(DEM2.z_masked))+300]);
        set(gca,'xlim',x_bounds,'ylim',y_bounds);
%         figure1b = figure; set(gcf,'position',[850 50 800 800]);
%         imagesc(IM.x,IM.y,IM.z); colormap gray;
%         set(gca,'ydir','normal'); set(gca,'xlim',x_bounds,'ylim',y_bounds);
        figure(figureA); 
        subplot(sub1); set(gca,'xlim',x_bounds,'ylim',y_bounds);
        subplot(sub2); set(gca,'xlim',x_bounds,'ylim',y_bounds);
        clear xv yv;
        disp('Draw polygon around the overlapping portion of the fjord with open water');
        figure(figure1);
        [fjord_mask,xv,yv] = roipoly;
        fjord_mask = double(fjord_mask);
        x1 = nearestneighbour(xv',DEM1.x); y1 = nearestneighbour(yv',DEM1.y);
        fjord.dem1.x = DEM1.x(min(x1):max(x1)); fjord.dem1.y = DEM1.y(min(y1):max(y1));
        x2 = nearestneighbour(xv',DEM2.x); y2 = nearestneighbour(yv',DEM2.y);
        fjord.dem2.x = DEM2.x(min(x2):max(x2)); fjord.dem2.y = DEM2.y(min(y2):max(y2));
        z2 = fjord_mask.*DEM2.z_masked; z2(z2 == 0) = NaN;
        fjord.dem2.z = z2(min(y2):max(y2),min(x2):max(x2));
        fjord.dem1.z = DEM1.z_masked(min(y1):max(y1),min(x1):max(x1));
        fjord_mask_crop = fjord_mask(min(y2):max(y2),min(x2):max(x2));
        clear xv yv;
        
        %interpolate to the same grid
        [x_grid_dem1,y_grid_dem1] = meshgrid(fjord.dem1.x,fjord.dem1.y);
        [x_grid_dem2,y_grid_dem2] = meshgrid(fjord.dem2.x,fjord.dem2.y);
        fjord.dem1.z_interp = interp2(x_grid_dem1,y_grid_dem1,fjord.dem1.z,x_grid_dem2,y_grid_dem2);
        fjord.dem1.z_interp = fjord_mask_crop.*fjord.dem1.z_interp; fjord.dem1.z_interp(fjord.dem1.z_interp == 0) = NaN;
        
        %plot the fjord data
        close all; drawnow;
        cmin = nanmedian(min(fjord.dem2.z));
        figure7 = figure;
        imagesc(fjord.dem2.x,fjord.dem2.y,fjord.dem1.z_interp);
        set(gca,'ydir','normal'); colorbar;
        set(gca,'clim',[cmin, cmin+60]);
        set(figure7,'position',[50 550 500 500]);
        figure8 = figure;
        imagesc(fjord.dem2.x,fjord.dem2.y,fjord.dem2.z);
        set(gca,'ydir','normal'); colorbar;
        set(gca,'clim',[cmin, cmin+60]);
        set(figure8,'position',[550 550 500 500]);
        
        %calculate statistics for the fjord offset
        fjord.dem_diff.map = fjord.dem1.z_interp - fjord.dem2.z;
        fjord.dem_diff.mean = nanmean(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        fjord.dem_diff.stdev = nanstd(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        fjord.dem_diff.median = nanmedian(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        fjord.dem_diff.mad = mad(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        
        %filter the fjord min elevation data
        fjord_filter = fjord.dem_diff.map >= fjord.dem_diff.median+1.48*fjord.dem_diff.mad | fjord.dem_diff.map <= fjord.dem_diff.median-1.48*fjord.dem_diff.mad;
        fjord.dem_diff.map(fjord_filter) = NaN;
        fjord.dem_diff.mean = nanmean(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        fjord.dem_diff.stdev = nanstd(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        fjord.dem_diff.median = nanmedian(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        fjord.dem_diff.mad = mad(fjord.dem_diff.map(~isnan(fjord.dem_diff.map)));
        disp(['Sea level offset: ',num2str(fjord.dem_diff.mean), 'm']);
        
        %save the fjord offset data
        save([dir_output,'fjord_offset_',DEM1.time,'-',DEM2.time,'.mat'],'fjord');
        
        % calculate the sea level adjustment using bedrock to adjust DEM
        % bias & the tidal model is desired
        if nargin == 8
            disp('Extract tidal change between acquisition dates');
            early_time = datenum(DEM1.time,'yyyymmddHHMMSS');
            late_time = datenum(DEM2.time,'yyyymmddHHMMSS');
            SD_time = [early_time:(late_time-early_time)/10:late_time];
            cd([dir_TMD,'DATA/']);
            [tidal_z,~] = tmd_tide_pred('Model_CATS2008',SD_time,fjord_lat,fjord_lon,'z'); %change in tide at the fjord mouth
            DEM1.date = SD_time(1); DEM2.date = SD_time(end);
            fjord.tidal_change = tidal_z(1)-tidal_z(end); %tidal change between image acquisition times
            if isnan(fjord.tidal_change)
                fjord.tidal_change = 0;
            end
            
            %coregister
            mean_bedrock_offset = 0; %assume bedrock offset is zero
            fjord.dem_diff.br_tide_adjust = mean_bedrock_offset + fjord.tidal_change;
            DEM2.z_masked_br_tide_adjust = DEM2.z_masked + fjord.dem_diff.br_tide_adjust;
            disp(['Sea level offset (bedrock & tide adjusted): ',num2str(fjord.dem_diff.br_tide_adjust), 'm']);
        end
        
        %coregister using the fjord minimum elevations
        close all;
        fjord.dem_diff.sl_adjust = fjord.dem_diff.median;
        DEM2.z_masked_sl_adjust = DEM2.z_masked + fjord.dem_diff.sl_adjust;
        disp(['Sea level offset (sea level adjusted): ',num2str(fjord.dem_diff.sl_adjust), 'm']);
        save([dir_output,'fjord_offset_',DEM1.time,'-',DEM2.time,'.mat'],'fjord');
        
    end
    
    
    %save the DEMs
    disp('Resaving DEMs with the sea level adjustment info added to their structures');
    DEM=DEM1; save([dir_output,DEM1.filename,'.mat'],'DEM','-v7.3','-append');
    clear DEM; DEM=DEM2; save([dir_output,DEM2.filename,'.mat'],'DEM','-v7.3','-append');
    clear DEM
    
end

end