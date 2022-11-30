function [DEM,IM] = convert_PGC_tifs_to_matfiles(DEM,dir_DEM,dir_output)
% Function to identify ortho image tifs corresponding to the DEMs produced
% by quality_check_DEMs.m
% Ellyn Enderlin (ellynenderlin@boisestate.edu)
% Fall 2022
% 
% INPUTS:   DEM             structure variable containing DEM filename 
%                           ([region_abbrev]-[image_capture_time]) 
%                           and image capture time (YYYYMMDD)
%           dir_DEM         directory where ortho .tif files exist
%           dir_output      directory where all output files will be placed
%
% OUTPUTS:  IM              8-bit structure variable of image file

%display a few outputs to check code
disp(['DEM path = ',dir_DEM]);
disp(['DEM date = ',DEM.filename(end-11:end-4)]);

%find the geotiffs for the specified date
datefiles = dir([dir_DEM,'SETSM*',DEM.filename(end-11:end-4),'*.tif']);
disp(['Number of SETSM files for that date: ',num2str(length(datefiles))]);

%load the DEM
dateref = 0;
for j = 1:length(datefiles)
    if contains(datefiles(j).name,'dem.tif')
        disp('loading DEM geotiff');
        [A,R] = readgeoraster([dir_DEM,datefiles(j).name]);
        if dateref == 0 %load the first or only file for this date
            DEM.x = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
            DEM.y = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;
            DEM.z = single(A); DEM.z(DEM.z == -9999) = NaN;
            clear A R;
        else
            %create grids for the DEM mosaic
            [xgrid1,ygrid1] = meshgrid(DEM.x,DEM.y); 
            [xgrid2,ygrid2] = meshgrid(R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX,...
                R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY); 
            demx_lims = [min([DEM.x(1),R.XWorldLimits(1)+0.5*R.CellExtentInWorldX]) max([DEM.x(2),R.XWorldLimits(2)-0.5*R.CellExtentInWorldX])];
            demy_lims = [min([DEM.y(1),R.YWorldLimits(2)-0.5*R.CellExtentInWorldY]) max([DEM.y(2),R.YWorldLimits(1)+0.5*R.CellExtentInWorldY])];
            [xgrid_mosaic,ygrid_mosaic] = meshgrid(demx_lims(1):R.CellExtentInWorldX:demx_lims(2),demy_lims(1):-R.CellExtentInWorldY:demy_lims(2));
            
            %create the mosaic
            dem_mosaic = NaN(size(xgrid));
            dem1 = interp2(xgrid1,ygrid1,double(DEM.z),xgrid_mosaic,ygrid_mosaic);
            dem2 = interp2(xgrid2,ygrid2,double(A),xgrid_mosaic,ygrid_mosaic);
            dem_mosaic = mean(cat(3,dem_mosaic,dem1,dem2),3,'omitnan');
            
            %replace original DEM with mosaic
            clear DEM; 
            DEM.x = demx_lims(1):R.CellExtentInWorldX:demx_lims(2);
            DEM.y = demy_lims(1):-R.CellExtentInWorldY:demy_lims(2);
            DEM.z = single(dem_mosaic);
            clear A R xgrid* ygrid* demx_lims demy_lims dem_mosaic dem1 dem2;
        end
        dateref = 1; %flag that there is more than one file for this date
    end
end



%create a mask
data_mask = zeros(size(DEM.z));
for i = 1:size(DEM.z,1)
    startval = find(~isnan(DEM.z(i,:))==1,1,'first');
    endval = find(~isnan(DEM.z(i,:))==1,1,'last');
    if ~isempty(startval)
        data_mask(i,startval:endval) = 1;
    end
    clear startval endval;
end
DEM.mask = data_mask; %removed Y.z_masked

%save DEM as a mat-file
save([dir_output,DEM.filename,'.mat'],'DEM','-v7.3'); % re-save the DEM data
disp('DEM saved');

%load the 8-bit image and save as a mat-file
dateref = 0;
for j = 1:length(datefiles)
    if contains(datefiles(j).name,'ortho.tif')
        disp('loading image geotiff');
        [A,R] = readgeoraster([dir_DEM,datefiles(j).name]);
        if dateref == 0 %load the first or only file for this date
            %create image coordinate vectors
            IM.x = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
            IM.y = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;
            
            %normalize the image
            A(A==-9999) = NaN;
            B = single(A); B(B == -32768) = NaN; clear A;
            C = B-min(B(~isnan(B))); clear B;
            D = C./max(C(~isnan(C))); clear C;
            IM.z = 255*D; 

            clear D R;
        else
            %normalize the image
            A(A==-9999) = NaN;
            B = single(A); B(B == -32768) = NaN; clear A;
            C = B-min(B(~isnan(B))); clear B;
            D = C./max(C(~isnan(C))); clear C;
            
            %create grids for the image mosaic
            [xgrid1,ygrid1] = meshgrid(IM.x,IM.y);
            [xgrid2,ygrid2] = meshgrid(R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX,...
                R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY);
            imx_lims = [min([IM.x(1),R.XWorldLimits(1)+0.5*R.CellExtentInWorldX]) max([IM.x(2),R.XWorldLimits(2)-0.5*R.CellExtentInWorldX])];
            imy_lims = [min([IM.y(1),R.YWorldLimits(2)-0.5*R.CellExtentInWorldY]) max([IM.y(2),R.YWorldLimits(1)+0.5*R.CellExtentInWorldY])];
            [xgrid_mosaic,ygrid_mosaic] = meshgrid(imx_lims(1):R.CellExtentInWorldX:imx_lims(2),imy_lims(1):-R.CellExtentInWorldY:imy_lims(2));
            
            %create the mosaic
            im_mosaic = NaN(size(xgrid));
            im1 = interp2(xgrid1,ygrid1,double(IM.z),xgrid_mosaic,ygrid_mosaic);
            im2 = interp2(xgrid2,ygrid2,double(D),xgrid_mosaic,ygrid_mosaic);
            im_mosaic = mean(cat(3,im_mosaic,im1,im2),3,'omitnan');
            
            %replace original DEM with mosaic
            clear IM;
            IM.x = imx_lims(1):R.CellExtentInWorldX:imx_lims(2);
            IM.y = imy_lims(1):-R.CellExtentInWorldY:imy_lims(2);
            IM.z = single(im_mosaic);
            clear D R xgrid* ygrid* imx_lims imy_lims im_mosaic im1 im2;
        end
        dateref = 1; %flag that there is more than one file for this date
    end
end


%create a mask for the image then infill holes inside data domain
im_mask = zeros(size(IM.z));
for i = 1:size(IM.z,1)
    startval = find(~isnan(IM.z(i,:))==1,1,'first');
    endval = find(~isnan(IM.z(i,:))==1,1,'last');
    if ~isempty(startval)
        im_mask(i,startval:endval) = 1;
    end
    clear startval endval;
end

%save the image
IM.mask = im_mask; %removed IM.z_masked
save([dir_output,DEM.filename(1:end-3),'orthoimage.mat'],'IM','-v7.3');
disp('Image saved');

end