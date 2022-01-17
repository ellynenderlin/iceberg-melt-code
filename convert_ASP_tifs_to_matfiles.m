function [DEM,IM] = convert_ASP_tifs_to_matfiles(DEM,dir_DEM,dir_output)
% Function to read the ASP tifs & save as mat-files
% Ellyn Enderlin & Mariama Dryak 
% Slightly reformatted by Rainey Aberle, Fall 2021
% 
% INPUTS:   DEM             structure variable containing DEM filename 
%                           ([region_abbrev]_[image_capture_time]) 
%                           and image capture time (YYYYMMDDhhssmm)
%           dir_DEM         directory where DEM .tif files exist
%           dir_output      directory where all output files will be placed
%
% OUTPUTS:  DEM             DEM structure variable with new fields
%           IM              8-bit structure variable of image file

%load the DEM
[A,R] = readgeoraster([dir_DEM,DEM.filename,'-DEM.tif']);
DEM.x = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
DEM.y = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;
DEM.z = single(A); DEM.z(DEM.z == -9999) = NaN;

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
save([dir_output,DEM.filename,'-DEM.mat'],'DEM','-v7.3'); % re-save the DEM data
disp('DEM saved');
clear A R;

%load the 8-bit image and save as a mat-file
[A,R] = readgeoraster([dir_DEM,DEM.filename,'-orthoimage.tif']); A(A==-9999) = NaN;
IM.x = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
IM.y = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;

%normalize the image
B = single(A); B(B == -32768) = NaN; clear A;
C = B-min(B(~isnan(B))); clear B;
D = C./max(C(~isnan(C))); clear C;

%create a mask for the image then infill holes inside data domain
im_mask = zeros(size(D));
for i = 1:size(D,1)
    startval = find(~isnan(D(i,:))==1,1,'first');
    endval = find(~isnan(D(i,:))==1,1,'last');
    if ~isempty(startval)
        im_mask(i,startval:endval) = 1;
    end
    clear startval endval;
end

%save the image
IM.z = 255*D; IM.mask = im_mask; %removed IM.z_masked
save([dir_output,DEM.filename,'-orthoimage.mat'],'IM','-v7.3');
disp('Image saved');

end