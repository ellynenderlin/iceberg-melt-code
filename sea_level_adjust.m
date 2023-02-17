function [A,B,sl_offset1,sl_offset2,cmin_early,cmin_late,flag] = sea_level_adjust(DEM1_time,DEM2_time,IM1,IM2,dir_output,A,B,S,sl_early,sl_late,fjord_offset,fjordz_o,fjordz_f,iceberg_no,region_abbrev,cmax,DEMo_pos,DEMf_pos,imo_pos,imf_pos)
% Function to adjust elevations using the estimated sea level around icebergs
% Ellyn Enderlin & Rainey Aberle
% Last edit: 17 Feb. 2023
%
% INPUTS:   DEM1_time       DEM1 mean image capture time
%           DEM2_time       DEM2 mean image capture time
%           IM1             structure variable containing earlier
%                               orthoimage info
%           IM2             structure variable containing later orthoimage
%                               info
%           dir_output      directory where output files will be saved
%           A               DEM1 cropped to iceberg region 
%           B               DEM2 cropped to iceberg region
%           S               iceberg shapefile
%           sl_early        sea level in earlier DEM 
%           sl_late         sea level in later DEM
%           fjord_offset    median offset of fjord between DEMs
%           fjordz_o        fjord elevation adjustment for DEM1
%           fjordz_f        fjord elevation adjustment for DEM2
%           iceberg_no      iceberg number 
%           region_abbrev   region abbreviation used in image files
%           cmax            
%           DEMo_pos        position of figure for plotting DEM1
%           DEMf_pos        position of figure for plotting DEM1
%           imo_pos         position of figure for plotting IM1
%           imf_pos         position of figure for plotting IM2
%
% OUTPUTS:  A               DEM1 cropped to iceberg region
%           B               DEM2 cropped to iceberg region
%           sl_offset1      sea level offset for DEM1
%           sl_offset2      sea level offset for DEM2
%           cmin_early      median vertical bias in earlier image
%           cmin_late       median vertical bias in later image
%           flag            flag for adjustment values (0=good, 1=unreliable)
elev_cmap = cmocean('thermal',10001); elev_cmap(1,:) = [1 1 1]; 

disp('Use nearby ice-free regions to estimate vertical bias between DEMs. Avoid shadows!!!');
close all; drawnow;

zmins2 = min(A.z_elpsd_adjust); ymins2 = min(B.z_elpsd_adjust);
zmins2(zmins2<-100 | zmins2>100) = NaN; ymins2(ymins2<-100 | ymins2>100) = NaN;
zmins2(zmins2<(nanmedian(zmins2)-3*1.4826*mad(zmins2,1)) | zmins2>(nanmedian(zmins2)+3*1.4826*mad(zmins2,1))) = NaN;
ymins2(ymins2<(nanmedian(ymins2)-3*1.4826*mad(ymins2,1)) | ymins2>(nanmedian(ymins2)+3*1.4826*mad(ymins2,1))) = NaN;
cmin_early = nanmedian(zmins2); cmin_late = nanmedian(ymins2);
% cmax_early = max(max(A.z_elpsd_adjust)); cmax_late = max(max(B.z_elpsd_adjust));
if isinf(cmin_early); cmin_early = -3; end
if isinf(cmin_late); cmin_late = -3; end
% if isinf(cmax_early); cmax_early = 50; end
% if isinf(cmax_late); cmax_late = 50; end
if isnan(cmin_early); cmin_early = -3; end
if isnan(cmin_late); cmin_late = -3; end
% if isnan(cmax_early); cmax_early = 50; end
% if isnan(cmax_late); cmax_late = 50; end
% cmax = max([cmax_early;cmax_late]);

%show the bedrock & tide adjustment for comparison
% uniform_adjust = fjord_offset+nansum(cat(3,early_elev_offset,late_elev_offset),3);
if length(fjord_offset)==2
    disp(['Regional fjord elevation-based adjustment = ', num2str(fjord_offset(1)),'m']);
    disp(['Bedrock elevation- & tide model-based adjustment = ', num2str(fjord_offset(2)),'m']);
else
    disp(['Regional fjord elevation-based adjustment = ', num2str(fjord_offset(1)),'m']);
end

q=1;
while q %use ice-free polygons to remove DEM bias (pick twice)
    %use ice-free regions to adjust sea level
    figure1 = figure; set(figure1,'position',DEMo_pos); %set(figure1,'position',[50 800 800 700]);
    imagesc(A.x,A.y,A.z_elpsd_adjust-cmin_early); hold on; axis xy equal; set(gca,'clim',sort([0 cmax-cmin_early]),'fontsize',14);
    colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
    title(['Early date: ',num2str(DEM1_time(1:8))],'fontsize',16);
    xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
    plot(S.X,S.Y,'-*k','linewidth',2,'markersize',4); hold on; %iceberg ROI
    plot(sl_early.X,sl_early.Y,'--y','linewidth',2,'markersize',2); hold on; %ice-free guide ROI
    set(gca,'xlim',xlims,'ylim',ylims);
    figureA = figure; set(figureA,'position',imo_pos); %set(figureA,'position',[850 800 800 700]);
%     imagesc(im1.x,im1.y,im1.z_adjust); axis xy equal; colormap gray; hold on;
    imagesc(IM1.x,IM1.y,IM1.z); axis xy equal; colormap gray; hold on;
    plot(S.X,S.Y,'-*r','linewidth',2,'markersize',4); hold on;
    set(gca,'xlim',xlims,'ylim',ylims);
    title(['Early date: ',num2str(DEM1_time(1:8))],'fontsize',16);
    figure2 = figure; set(figure2,'position',DEMf_pos); %set(figure2,'position',[50 50 800 700]);
    imagesc(B.x,B.y,B.z_elpsd_adjust-cmin_late); hold on; axis xy equal; set(gca,'clim',sort([0 cmax-cmin_late]),'fontsize',14);
    colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
    title(['Late date: ',num2str(DEM2_time(1:8))],'fontsize',16);
    xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
    plot(sl_late.X,sl_late.Y,'--y','linewidth',2,'markersize',2); hold on; %ice-free guide ROI
    set(gca,'xlim',xlims,'ylim',ylims);
    figureB = figure; set(figureB,'position',imf_pos); %set(figureB,'position',[850 50 800 700]);
%     imagesc(im2.x,im2.y,im2.z_adjust); axis xy equal; colormap gray; hold on;
    imagesc(IM2.x,IM2.y,IM2.z); axis xy equal; colormap gray; hold on;
    set(gca,'xlim',xlims,'ylim',ylims);
    title(['Late date: ',num2str(DEM2_time(1:8))],'fontsize',16);
    figure(figureA); figure(figure1);
    disp(['Picking ice-free ROIs: Iteration #',num2str(q),' of 5']);
    
    %outline an ice-free region with reliable elevations near the iceberg
    figure(figure1);
    disp('Draw an ice-free ROI in the early DEM using the initial guess as a guide');
    [poly1,s1x,s1y] = roipoly;
    icefree_area1 = poly1.*A.z_elpsd_adjust; 
    icefree_area1(icefree_area1 == 0) = NaN; %set everything already at sea level = NaN
    figure(figureA);
    plot(s1x,s1y,'--b','linewidth',2); hold on;
    figure(figure2);
    disp('Draw an ice-free ROI in the late DEM using the initial guess as a guide');
    [poly2,s2x,s2y] = roipoly;
    icefree_area2 = poly2.*B.z_elpsd_adjust; 
    icefree_area2(icefree_area2 == 0) = NaN; %set everything already at sea level = NaN
    figure(figureB);
    plot(s2x,s2y,'--b','linewidth',2); hold on;
    
    %adjust the iceberg DEM using the median elevation of ice-free pixels
    sl_dz1 = icefree_area1(~isnan(icefree_area1)); sl_dz1(sl_dz1>nanmedian(sl_dz1)+3*(1.4826*mad(sl_dz1,1)) | sl_dz1<nanmedian(sl_dz1)-3*(1.4826*mad(sl_dz1,1))) = NaN;
    sl_dz2 = icefree_area2(~isnan(icefree_area2)); sl_dz2(sl_dz2>nanmedian(sl_dz2)+3*(1.4826*mad(sl_dz2,1)) | sl_dz2<nanmedian(sl_dz2)-3*(1.4826*mad(sl_dz2,1))) = NaN;
    sl_offset1 = nanmedian(sl_dz1); if isnan(sl_offset1); sl_offset1 = 0; end
    sl_offset2 = nanmedian(sl_dz2); if isnan(sl_offset2); sl_offset2 = 0; end
    A.z_local_adjust = A.z_elpsd_adjust - sl_offset1;
    disp('Re-plotting DEMs coregistered with local sea level pixels');
    clf(figure1); figure(figure1);
    imagesc(A.x,A.y,A.z_local_adjust); hold on; axis xy equal; set(gca,'clim',[0 (cmax-(sl_offset2 + sl_offset1)/2)],'fontsize',14);
    colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
    title(['Early date: ',num2str(DEM1_time(1:8))],'fontsize',16);
    plot(S.X,S.Y,'-*k','linewidth',2,'markersize',4); hold on; %iceberg ROI
    B.z_local_adjust = B.z_elpsd_adjust - sl_offset2;
    clf(figure2); figure(figure2);
    imagesc(B.x,B.y,B.z_local_adjust); hold on; axis xy equal; set(gca,'clim',[0 (cmax-(sl_offset2 + sl_offset1)/2)],'fontsize',14);
    colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
    title(['Late date: ',num2str(DEM2_time(1:8))],'fontsize',16);
    drawnow;
    
    %compare the early and late DEM adjustments
    local_adjustment = -sl_offset2 + sl_offset1;
    disp(['Sea level-based adjustment = ',num2str(local_adjustment),' m']);
    if length(fjord_offset)==2
        adjust_diff = local_adjustment - fjord_offset(1);
        disp(['Local sea level - regional fjord elevation = ', num2str(adjust_diff),'m']);
        adjust_diff = local_adjustment - fjord_offset(2);
        disp(['Local sea level - bedrock elevation- & tide model-based adjustment = ', num2str(adjust_diff),'m']);
    else
        adjust_diff = local_adjustment - fjord_offset(1);
        disp(['Local sea level - regional fjord elevation = ', num2str(adjust_diff),'m']);
    end
    
    %check the agreement between the two bias estimation techniques
    disp('If the adjustment is good:');
    disp('1) random noise will cause some pixels at sea level to be slightly <0 (white) & some slightly >0 (dark blue)');
    disp('2) the offset between regional & local sea level estimates will probably be small');
    prompt = 'Are the sea level corrections reasonable (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        %write shapefiles for the good ice-free polygons
        cd([dir_output,'/',DEM1_time,'-',DEM2_time,'/icefree_rois/']);
        s.Geometry = 'Polygon';
        s.BoundingBox = [min(s1x) min(s1y); max(s1x) max(s1y)];
        s.X = double(s1x'); s.Y = double(s1y');
        s.Name = ['iceberg',num2str(iceberg_no)];
        shapefile_name = [region_abbrev,'_',num2str(DEM1_time),'_icefree',num2str(iceberg_no)];
        shapewrite(s,shapefile_name);
        clear s;
        s.Geometry = 'Polygon';
        s.BoundingBox = [min(s2x) min(s2y); max(s2x) max(s2y)];
        s.X = double(s2x'); s.Y = double(s2y');
        s.Name = ['iceberg',num2str(iceberg_no)];
        shapefile_name = [region_abbrev,'_',num2str(DEM2_time),'_icefree',num2str(iceberg_no)];
        shapewrite(s,shapefile_name);
        flag = 0; %flag it as having good adjustment values

        clear zmins2 ymins2;
        zmins2 = min(A.z_local_adjust); ymins2 = min(B.z_local_adjust);
        zmins2(zmins2<-100 | zmins2>100) = NaN; ymins2(ymins2<-100 | ymins2>100) = NaN;
        zmins2(zmins2<(nanmedian(zmins2)-3*1.4826*mad(zmins2,1)) | zmins2>(nanmedian(zmins2)+3*1.4826*mad(zmins2,1))) = NaN;
        ymins2(ymins2<(nanmedian(ymins2)-3*1.4826*mad(ymins2,1)) | ymins2>(nanmedian(ymins2)+3*1.4826*mad(ymins2,1))) = NaN;
        cmin_early = floor(nanmedian(zmins2)); cmin_late = floor(nanmedian(ymins2));
        
        break
    else
        if q<5
            disp('Repick ice-free regions for sl-adjustment');
            q = q+1;
            close all; drawnow;
            
        elseif q == 5
            disp('Cannot get good elevation adjustment from ice-free regions, using regional fjord adjustment.');
            sl_offset1 = nanmedian(min(fjordz_o)); sl_offset2 = nanmedian(min(fjordz_f));
            
            %write shapefiles for the last attempt at ice-free polygons
            cd([dir_output,'/',DEM1_time,'-',DEM2_time,'/icefree_rois']);
            s.Geometry = 'Polygon';
            s.BoundingBox = [min(s1x) min(s1y); max(s1x) max(s1y)];
            s.X = double(s1x'); s.Y = double(s1y');
            s.Name = ['iceberg',num2str(iceberg_no)];
            shapefile_name = [region_abbrev,'_',num2str(DEM1_time),'_icefree',num2str(iceberg_no)];
            shapewrite(s,shapefile_name);
            clear s;
            s.Geometry = 'Polygon';
            s.BoundingBox = [min(s2x) min(s2y); max(s2x) max(s2y)];
            s.X = double(s2x'); s.Y = double(s2y');
            s.Name = ['iceberg',num2str(iceberg_no)];
            shapefile_name = [region_abbrev,'_',num2str(DEM2_time),'_icefree',num2str(iceberg_no)];
            shapewrite(s,shapefile_name);
            flag = 1; %note that the adjustment is not reliable
            
            %replot figures with regional elevation adjustments for coregistration
            clf(figure1); figure(figure1);
            A.z_local_adjust = A.z_elpsd_adjust - sl_offset1;
            imagesc(A.x,A.y,A.z_reg_adjust); hold on; axis xy equal; set(gca,'clim',[0 (cmax-(sl_offset2 + sl_offset1)/2)],'fontsize',14);
            colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
            title(['Early date: ',num2str(DEM1_time(1:8))],'fontsize',16);
            plot(S.X,S.Y,'-*k','linewidth',2,'markersize',4); hold on; %iceberg ROI
            B.z_local_adjust = B.z_elpsd_adjust - sl_offset2;
            clf(figure2); figure(figure2);
            imagesc(B.x,B.y,B.z_reg_adjust); hold on; axis xy equal; set(gca,'clim',[0 (cmax-(sl_offset2 + sl_offset1)/2)],'fontsize',14);
            colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
            title(['Late date: ',num2str(DEM2_time(1:8))],'fontsize',16);
            
            zmins2 = min(A.z_local_adjust); ymins2 = min(B.z_local_adjust);
            zmins2(zmins2<-100 | zmins2>100) = NaN; ymins2(ymins2<-100 | ymins2>100) = NaN;
            zmins2(zmins2<(nanmedian(zmins2)-3*1.4826*mad(zmins2,1)) | zmins2>(nanmedian(zmins2)+3*1.4826*mad(zmins2,1))) = NaN;
            ymins2(ymins2<(nanmedian(ymins2)-3*1.4826*mad(ymins2,1)) | ymins2>(nanmedian(ymins2)+3*1.4826*mad(ymins2,1))) = NaN;
            cmin_early = floor(nanmedian(zmins2)); cmin_late = floor(nanmedian(ymins2));
            break
        end
    end
    
end

end
