function [SL] = update_or_remove_Antarctic_iceberg_meltrates(DEM1,DEM2,IM1,IM2,region_name,region_abbrev,bad_bergs,dir_output,dir_iceberg,dir_code)
% Function to convert iceberg elevation change to melt rates in Antarctica
% Ellyn Enderlin & Mariama Dryak
% Slightly Modified by Rainey Aberle, Fall 2021
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
% OUTPUTS:  SL structure containing iceberg melt volume & melt rate data
close all; drawnow;

%specify densities
rho_sw = 1026;  rho_sw_err = 2; %kg m^-3
rho_i = 917; rho_i_err = 10; %kg m^-3

%load the saved iceberg data
DEM1_date = DEM1.time(1:8); DEM2_date = DEM2.time(1:8);
load([dir_iceberg,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat']);
berg_no = [];
for i = 1:length(SL)
    berg_name(i,:) = SL(i).name;
    berg_no = [berg_no; SL(i).name(end-1:end)];
end

%load the firn density info
cd(dir_output);
load(['firn_data/',region_name,'_density_data.mat']);
density_levels = [700 750 800 830]; density_depths = [density.sevhun density.sevfif density.eighthun density.eightthir];
density_levels = [density_levels 917]; density_depths = [density_depths density.nineseventeen];
cd(dir_iceberg);

%only repeat calculations for troublesome icebergs
berg_refs = [];
if ~isempty(bad_bergs)
    for j = 1:length(bad_bergs)
        if length(num2str(bad_bergs(j)))
            bad_ref = matches(string(berg_no),['0',num2str(bad_bergs(j))]);
        else
            bad_ref = matches(string(berg_no),num2str(bad_bergs(j)));
        end
        berg_refs = [berg_refs; bad_ref];
        clear bad_ref;
    end
end

%calculate the time separation between DEMs in terms of 
%decimal years (ddays) & decimal days (days)
to = DEM1.time; tf = DEM2.time;
if mod(str2num(to(1:4)),4)==0; doyso=366; modayso = [31 29 31 30 31 30 31 31 30 31 30 31]; else doyso=365; modayso = [31 28 31 30 31 30 31 31 30 31 30 31]; end
if mod(str2num(tf(1:4)),4)==0; doysf=366; modaysf = [31 29 31 30 31 30 31 31 30 31 30 31]; else doysf=365; modaysf = [31 28 31 30 31 30 31 31 30 31 30 31]; end
doyo = sum(modayso(1:str2num(to(5:6))))-31+str2num(to(7:8)); doyf = sum(modaysf(1:str2num(tf(5:6))))-31+str2num(tf(7:8));
if str2num(tf(1:4)) == str2num(to(1:4))
    ddays = doyf-doyo+1;
elseif str2num(tf(1:4)) - str2num(to(1:4)) == 1
    ddays = doyf + (doyso-doyo)+1;
else
    years = str2num(to(1:4)):1:str2num(tf(1:4));
    for k = 1:length(years)
        if mod(years(k),4)==0
            doys(k)=366;
        else
            doys(k) = 365;
        end
    end
    
        
    %calculate the sum of days differently if during a leap year
    if doyo > sum(modayso(1:2))
        ddays = doyf + sum(doys(2:end-1)) + (365-doyo)+1;
    else
        ddays = doyf + sum(doys(2:end-1)) + (366-doyo)+1;
    end
end
hrs_o = ((str2num(to(13:14))/(60*60*24))+(str2num(to(11:12))/(60*24))+(str2num(to(9:10))/24));
hrs_f = ((str2num(tf(13:14))/(60*60*24))+(str2num(tf(11:12))/(60*24))+(str2num(tf(9:10))/24));
dhrs = hrs_f - hrs_o;
dt = ddays + dhrs;

%extract DEM pixel areas
DEM1_pixel_area = abs(DEM1.x(1)-DEM1.x(2)).*abs(DEM1.y(1)-DEM1.y(2)); DEM2_pixel_area = abs(DEM2.x(1)-DEM2.x(2)).*abs(DEM2.y(1)-DEM2.y(2)); %square meters

%set image bounds & crop
%EARLIER DATE
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
clear x1 x2 y1 y2 xlims ylims xmin xmax ymin ymax;
%LATER DATE
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
B.x = IM2.x(1,[xlims(1):xlims(2)]); B.y = IM2.y(1,[ylims(1):ylims(2)]);
B.z = double(IM2.z(ylims(1):ylims(2),xlims(1):xlims(2)));
if mean(B.z(~isnan(B.z)),'omitnan') < (1/3)*255 %image is dark
    z_adjust = imadjust(B.z./255,[],[],0.25);
    B.z = z_adjust; clear z_adjust;
end

%extract image pixel areas
im1_pixel_area = abs(A.x(1)-A.x(2)).*abs(A.y(1)-A.y(2)); im2_pixel_area = abs(B.x(1)-B.x(2)).*abs(B.y(1)-B.y(2)); %square meters
clear IM*;

%density data
%create density profiles
clear FAC;
FAC(1) = firnair.median; FAC(2) = firnair.median-firnair.uncert; FAC(3) = firnair.median+firnair.uncert; %estimate firn air content
ft = fittype('917-(917-a)*exp(-x/b)');[f,~] = fit(density_depths',density_levels',ft,'StartPoint',[350 firnair.median]); ci = confint(f); %create density profile
density_z = [0:1:1000];
density_profile = rho_i-(rho_i-f.a)*exp(-density_z/f.b); mindensity_profile = rho_i-(rho_i-ci(1,1))*exp(-density_z/ci(1,2)); maxdensity_profile = rho_i-(rho_i-ci(2,1))*exp(-density_z/ci(2,2));
%calculate wet density profile by flipping the shape of the exponential curve & compressing the range from (830-0) to (rho_sw-830)
wetdensity_profile = 830+((830-density_profile)./830).*(rho_sw-830); wetdensity_profile(ceil(density.eightthir)+1:end) = density_profile(ceil(density.eightthir)+1:end);
minwetdensity_profile = 830+((830-mindensity_profile)./830).*(rho_sw-830); minwetdensity_profile(ceil(density.eightthir)+1:end) = mindensity_profile(ceil(density.eightthir)+1:end);
maxwetdensity_profile = 830+((830-maxdensity_profile)./830).*(rho_sw-830); maxwetdensity_profile(ceil(density.eightthir)+1:end) = maxdensity_profile(ceil(density.eightthir)+1:end);


%update iceberg data for the icebergs for which you re-ran volume change calculations
cd(iceberg_dir);
for i = berg_refs
    if i < 10
        berg_number = [num2str(0),num2str(i)];
        load([iceberg',num2str(0),num2str(i),'_dz.mat']);
    else
        berg_number = num2str(i);
        load([iceberg',num2str(i),'_dz.mat']);
    end
    berg_no = strmatch([region_name,'iceberg',berg_number],berg_name);
    
    %incorporate data into a structure
    SL(berg_no).name = [region_name,'iceberg',num2str(berg_number)];
    xo = []; yo = []; xf = []; yf = [];
    for j = 1:10
        xo = cat(1,xo,IB(j).vertices.xo); yo = cat(1,yo,IB(j).vertices.yo);
        xf = cat(1,xf,IB(j).vertices.xf); yf = cat(1,yf,IB(j).vertices.yf);
    end
    SL(berg_no).initial.x = nanmean(xo,1); SL(berg_no).initial.y = nanmean(yo,1);
    SL(berg_no).final.x = nanmean(xf,1); SL(berg_no).final.y = nanmean(yf,1);
    poly1 = roipoly(DEM1.x,DEM1.y,DEM1.z_elpsd_adjust,SL(berg_no).initial.x,SL(berg_no).initial.y);
    mask1 = logical(poly1);
    SL(berg_no).initial.z = single(mask1.*DEM1.z_elpsd_adjust); SL(berg_no).initial.z(SL(berg_no).initial.z == 0) = NaN; SL(berg_no).initial.z(SL(berg_no).initial.z<0) = NaN;
    poly2 = roipoly(DEM2.x,DEM2.y,DEM2.z_elpsd_adjust,SL(berg_no).final.x,SL(berg_no).final.y);
    mask2 = logical(poly2);
    SL(berg_no).final.z = single(mask2.*DEM2.z_elpsd_adjust); SL(berg_no).final.z(SL(berg_no).final.z == 0) = NaN; SL(berg_no).final.z(SL(berg_no).final.z<0) = NaN;
    SL(berg_no).initial.coreg_z = IB(1).local_adjust_o; SL(berg_no).final.coreg_z = IB(1).local_adjust_f;
    if isfield(IB,'flag')
        SL(berg_no).flag = IB(1).flag;
    else
        SL(berg_no).flag = NaN;
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
    z = DEM1.z_elpsd_adjust+SL(berg_no).initial.coreg_z*ones(size(DEM1.z_elpsd_adjust));
    idx = nearestneighbour([A.x(1) A.x(end)],DEM1.x); idy = nearestneighbour([A.y(1) A.y(end)],DEM1.y);
    DEM_x = DEM1.x(idx(1):idx(2)); DEM_y = DEM1.y(idy(1):idy(2));
    DEM_z = z(idy(1):idy(2),idx(1):idx(2));
    imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    
    %measure the iceberg in the early image & DEM
    disp('Measure the iceberg aerial extent by drawing a polygon around the iceberg edge');
    cd_to_glacier_dir = ['cd ',glacier_dir]; eval(cd_to_glacier_dir);
    cd iceberg_data/iceberg_shapes
    disp(['Zooming-in & plotting the DEM ROI in the image and DEM for iceberg #',num2str(i)]);
    figure(figure1);
    %     imagesc(A.x,A.y,A.z); set(gca,'ydir','normal'); hold on;
    colormap gray; set(gca,'clim',[1.05*min(A.z(~isnan(A.z))) (0.95)*max(max(A.z))]);
    vxi = nearestneighbour(SL(berg_no).initial.x,A.x); vyi = nearestneighbour(SL(berg_no).initial.y,A.y);
    set(gca,'xlim',[min(A.x(vxi))-150 max(A.x(vxi))+150],'ylim',[min(A.y(vyi))-150 max(A.y(vyi))+150]);
    plot(A.x(vxi),A.y(vyi),'--r','linewidth',2);
    figure(figure2);
    %     imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    set(gca,'xlim',[min(SL(berg_no).initial.x)-150 max(SL(berg_no).initial.x)+150],'ylim',[min(SL(berg_no).initial.y)-150 max(SL(berg_no).initial.y)+150]);
    plot(SL(berg_no).initial.x,SL(berg_no).initial.y,'-k','linewidth',2); hold on; plot(SL(berg_no).initial.x,SL(berg_no).initial.y,'--w','linewidth',2); hold on;
    
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
    SL(berg_no).initial.x_perim = x; SL(berg_no).initial.y_perim = y;
    
    %modify the mask using the DEM if necessary
    disp('Now modify the iceberg mask needed to extract elevations from the entire iceberg area');
    clear DEM_z_masked iceberg_DEMmask x y;
    figure(figure2);
    [iceberg_DEMmask,x,y] = roipoly;
    SL(berg_no).initial.x_perim = []; SL(berg_no).initial.y_perim = [];
    SL(berg_no).initial.x_perim = x; SL(berg_no).initial.y_perim = y;
    DEM_z_masked = iceberg_DEMmask.*DEM_z;
    DEM_z_masked(iceberg_DEMmask == 0) = 0;
    clf(figure2); drawnow;
    imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    set(gca,'xlim',[min(SL(berg_no).initial.x)-150 max(SL(berg_no).initial.x)+150],'ylim',[min(SL(berg_no).initial.y)-150 max(SL(berg_no).initial.y)+150]);
    plot(SL(berg_no).initial.x,SL(berg_no).initial.y,'-k','linewidth',2); hold on; plot(SL(berg_no).initial.x,SL(berg_no).initial.y,'--w','linewidth',2); hold on;
    plot(SL(berg_no).initial.x_perim,SL(berg_no).initial.y_perim,'-w','linewidth',2); hold on; plot(SL(berg_no).initial.x_perim,SL(berg_no).initial.y_perim,'--m','linewidth',2); hold on;
    
    %use the mask to extract size & shape info
    disp('Pulling plan-view info');
    for k = 1:size(x,1)-1
        dist(k) = ((x(k)-x(k+1))^2+(y(k)-y(k+1))^2)^(1/2);
    end
    dist(end+1) = ((x(end)-x(1))^2+(y(end)-y(1))^2)^(1/2);
    shapes_o = regionprops(iceberg_IMmask,'Area','MinorAxisLength','MajorAxisLength'); %use the image-derived mask to estimate area
    lmin_o = abs(A.x(2)-A.x(1))*shapes_o.MinorAxisLength;
    lmax_o = abs(A.x(2)-A.x(1))*shapes_o.MajorAxisLength;
    SL(berg_no).initial.width_min = lmin_o; SL(berg_no).initial.width_max = lmax_o;
    base_area = shapes_o.Area.*im1_pixel_area; SL(berg_no).initial.SA = base_area;
    SL(berg_no).initial.z_map = single(DEM_z_masked); SL(berg_no).initial.z_map(SL(berg_no).initial.z_map == 0) = NaN; SL(berg_no).initial.z_map(SL(berg_no).initial.z_map<0) = NaN;
    zo = SL(berg_no).initial.z_map(~isnan(SL(berg_no).initial.z_map));
    SL(berg_no).initial.z_median = nanmedian(zo); SL(berg_no).initial.z_mad = mad(zo,1);
    perimeter = sum(dist);
    SL(berg_no).initial.radius = perimeter/(2*pi);
    clear iceberg_IMmask;
    
    %iteratively estimate bulk density
        if SL(i).orientation == 0 %berg is a fragment or overturned so assume the full firn is saturated
        for k = 1:3 %k=1 is the best guess, k=[2,3] constraints uncertainties
            j=1;
            Hberg(j) = (rho_sw/(rho_sw-rho_i))*SL(i).initial.z_median;
            Hberg_err(j) = abs(Hberg(j)).*sqrt((rho_sw_err/rho_sw)^2 + ((1.4826*SL(i).initial.z_mad)/SL(i).initial.z_median)^2 + ((rho_sw_err^2 + rho_i_err^2)/(rho_sw-rho_i)^2));
            wet_ref = find(density_z<=density.eightthir,1,'last');
            if Hberg(j) > density_z(wet_ref)
                rho_prof = [wetdensity_profile(1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [minwetdensity_profile(1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxwetdensity_profile(1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
            else
                rho_prof = [wetdensity_profile(1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [minwetdensity_profile(1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxwetdensity_profile(1:ceil(Hberg(j))+1)];
            end
            rho_f(j) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
            rho_f_err(j) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]));
            clear rho_prof*;
            
            while j
                Hberg(j+1) = (rho_sw/(rho_sw-rho_f(j)))*SL(i).initial.z_median;
                Hberg_err(j+1) = abs(Hberg(j+1)).*sqrt(((abs(rho_sw/(rho_sw-rho_f(j)))*sqrt((rho_sw_err/rho_sw)^2 + (sqrt(rho_sw_err^2+rho_f_err(j)^2)/(rho_sw-rho_f(j)))^2))/(rho_sw/(rho_sw-rho_f(j))))^2 + ((1.4826*SL(i).initial.z_mad)/SL(i).initial.z_median)^2);
                if Hberg(j+1) > density_z(wet_ref)
                    rho_prof = [wetdensity_profile(1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [minwetdensity_profile(1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxwetdensity_profile(1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                else
                    rho_prof = [wetdensity_profile(1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [minwetdensity_profile(1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxwetdensity_profile(1:ceil(Hberg(j+1))+1)];
                end
                rho_f(j+1) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
                rho_f_err(j+1) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]));
                clear rho_prof*;

                if abs(rho_f(j+1)-rho_f(j)) < 0.25*rho_f_err(j+1)
                    if k == 1
                        SL(i).density_type = 'wet';
                        SL(i).initial.density = rho_f(j+1);
                        %propagated density uncertainty
                        SL(i).initial.density_uncert = rho_f_err(j);
                    else
                        %density uncertainty using FAC uncertainty bounds to iterate thickness & rho_f estimates
                        SL(i).initial_range.density(k-1) = rho_f(j+1);
                    end
                    clear Hberg rho_f* dry_ref wet_ref;
                    break
                else
                    j = j+1;
                end
%                 clear Hice Vice;
            end
        end
    else %berg is upright so only wet the firn below the waterline
        for k = 1:3 %k=1 is the best guess, k=[2,3] constraints uncertainties
            j=1;
            Hberg(j) = (rho_sw/(rho_sw-rho_i))*SL(i).initial.z_median;
            Hberg_err(j) = abs(Hberg(j)).*sqrt((rho_sw_err/rho_sw)^2 + ((1.4826*SL(i).initial.z_mad)/SL(i).initial.z_median)^2 + ((rho_sw_err^2 + rho_i_err^2)/(rho_sw-rho_i)^2));
            dry_ref = find(density_z<=SL(i).initial.z_median,1,'last'); wet_ref = find(density_z<=density.eightthir,1,'last');
            if Hberg(j) > density_z(wet_ref)
                rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
            else
                rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:ceil(Hberg(j))+1)];
            end
            rho_f(j) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
            rho_f_err(j) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]));
            clear rho_prof*;
            
            while j
                Hberg(j+1) = (rho_sw/(rho_sw-rho_f(j)))*SL(i).initial.z_median;
                Hberg_err(j+1) = abs(Hberg(j+1)).*sqrt(((abs(rho_sw/(rho_sw-rho_f(j)))*sqrt((rho_sw_err/rho_sw)^2 + (sqrt(rho_sw_err^2+rho_f_err(j)^2)/(rho_sw-rho_f(j)))^2))/(rho_sw/(rho_sw-rho_f(j))))^2 + ((1.4826*SL(i).initial.z_mad)/SL(i).initial.z_median)^2);
                if Hberg(j+1) > density_z(wet_ref)
                    rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                else
                    rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:ceil(Hberg(j+1))+1)];
                end
                rho_f(j+1) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
                rho_f_err(j+1) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]));
                clear rho_prof*;
                
                if abs(rho_f(j+1)-rho_f(j)) < 0.25*rho_f_err(j+1)
                    if k == 1
                        SL(i).density_type = 'dry';
                        SL(i).initial.density = rho_f(j+1);
                        %propagated density uncertainty
                        SL(i).initial.density_uncert = rho_f_err(j);
                    else
                        %density uncertainty using FAC uncertainty bounds to iterate thickness & rho_f estimates
                        SL(i).initial_range.density(k-1) = rho_f(j+1);
                    end
                    clear Hberg rho_f* dry_ref wet_ref;
                    break
                else
                    j = j+1;
                end
%                 clear Hice Vice;
            end
        end
    end
    
    %estimate the iceberg depth & submerged area using densities from the Ligtenberg FDM
    disp('Estimating thickness & volume');
    if SL(berg_no).orientation == 0 %overturned/fragment
        rho_f = sort([rho_i SL(berg_no).initial.density]); %figure out if bubble-free ice or water-saturated firn has a lower density
        for l = 1:2
            draft(l) = (rho_f(l)/(rho_sw-rho_f(l)))*SL(berg_no).initial.z_median;
            Hberg(l) = (rho_sw/(rho_sw-rho_f(l)))*SL(berg_no).initial.z_median;
        end
    else
        rho_f = sort(SL(berg_no).initial_range.density); %assume density constrained by model FAC range
        for l = 1:2
            draft(l) = (rho_f(l)/(rho_sw-rho_f(l)))*SL(berg_no).initial.z_median;
            Hberg(l) = (rho_sw/(rho_sw-rho_f(l)))*SL(berg_no).initial.z_median;
        end
    end
    SL(berg_no).initial_range.draft = sort(draft);
    lat_area = sort(draft).*sum(dist); SL(berg_no).initial_range.LA = lat_area;
    area = base_area*ones(size(lat_area)) + lat_area; SL(berg_no).initial_range.TA = area;
    SL(berg_no).initial.V = DEM1_pixel_area*(rho_sw/(rho_sw-SL(berg_no).initial.density))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(berg_no).initial.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(berg_no).initial_range.V(1) = DEM1_pixel_area*(rho_sw/(rho_sw-min(rho_f)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(berg_no).initial.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(berg_no).initial_range.V(2) = DEM1_pixel_area*(rho_sw/(rho_sw-max(rho_f)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(berg_no).initial.z_median*sum(sum(isnan(DEM_z_masked))));
    
    
    %save the iceberg outline
    disp('Writing shapefile containing early iceberg polygon information');
    S.Geometry = 'Polygon';
    S.BoundingBox = [min(x) min(y); max(x) max(y)];
    S.X = double(x'); S.Y = double(y');
    S.Name = ['iceberg',num2str(berg_number)];
    shapefile_name = ['WV_',num2str(to),'_icebergshape',num2str(berg_number)];
    cd([dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);
    shapewrite(S,shapefile_name);
    copyfile([dir_code,'general/antarctic_PSprojection.prj'],[dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/',shapefile_name,'.prj']);
    cd([dir_output,DEM1.time,'-',DEM2.time,'/']);
    clear S;
    
    %save the data
    disp('Saving data');
    save([dir_output,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat'],'SL','-v7.3');
    cd([dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);
    clear DEM_z_masked z idx idy x y iceberg_IMmask iceberg_DEMmask rho_f draft Hberg lat_area base_area dist area_mask masked_pixels image_xgrid image_ygrid shapes lmin* lmax* perimeter;
    clear x1 x2 y1 y2 xlims ylims xmin xmax ymin ymax DEM_x DEM_y DEM_z;
    disp('Advancing to the later date');
    close(figure1); close(figure2); drawnow;
    
    %plot the later DEM & image
    disp('Plotting the later image-DEM pair for the fjord');
    figure1 = figure;
    imagesc(B.x,B.y,B.z); hold on;
    colormap gray; set(gca,'clim',[1.05*min(B.z(~isnan(B.z))) (0.95)*max(max(B.z))]);
    set(gca,'ydir','normal'); hold on;
    set(gcf,'position',[50 100 600 600]);
    figure2 = figure;
    set(gcf,'position',[650 100 600 600]);
    z = DEM2.z_elpsd_adjust+SL(berg_no).final.coreg_z*ones(size(DEM2.z_elpsd_adjust));
    idx = nearestneighbour([B.x(1) B.x(end)],DEM2.x); idy = nearestneighbour([B.y(1) B.y(end)],DEM2.y);
    DEM_x = DEM2.x(idx(1):idx(2)); DEM_y = DEM2.y(idy(1):idy(2));
    DEM_z = z(idy(1):idy(2),idx(1):idx(2));
    imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    
    %measure the iceberg in the later image & DEM
    cd_to_glacier_dir = ['cd ',glacier_dir]; eval(cd_to_glacier_dir);
    cd iceberg_data/iceberg_shapes
    disp(['Zooming-in & plotting the DEM ROI in the image and DEM for iceberg #',num2str(i)]);
    figure(figure1);
    colormap gray; set(gca,'clim',[1.05*min(B.z(~isnan(B.z))) (0.95)*max(max(B.z))]);
    vxf = nearestneighbour(SL(berg_no).final.x,B.x); vyf = nearestneighbour(SL(berg_no).final.y,B.y);
    set(gca,'xlim',[min(SL(berg_no).final.x)-150 max(SL(berg_no).final.x)+150],'ylim',[min(SL(berg_no).final.y)-150 max(SL(berg_no).final.y)+150]);
    plot(B.x(vxf),B.y(vyf),'--r','linewidth',2);
    figure(figure2);
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    set(gca,'xlim',[min(SL(berg_no).final.x)-150 max(SL(berg_no).final.x)+150],'ylim',[min(SL(berg_no).final.y)-150 max(SL(berg_no).final.y)+150]);
    plot(SL(berg_no).final.x,SL(berg_no).final.y,'-k','linewidth',2); hold on; plot(SL(berg_no).final.x,SL(berg_no).final.y,'--w','linewidth',2); hold on;
    
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
    SL(berg_no).final.x_perim = x; SL(berg_no).final.y_perim = y;
    
    %modify the mask using the DEM if necessary
    disp('Now modify the iceberg mask need to extract elevations from the entire iceberg area');
    clear DEM_z_masked iceberg_DEMmask x y;
    figure(figure2);
    [iceberg_DEMmask,x,y] = roipoly;
    SL(berg_no).final.x_perim = []; SL(berg_no).final.y_perim = [];
    SL(berg_no).final.x_perim = x; SL(berg_no).final.y_perim = y;
    DEM_z_masked = iceberg_DEMmask.*DEM_z;
    DEM_z_masked(iceberg_DEMmask == 0) = 0;
    clf(figure2); drawnow;
    imagesc(DEM_x,DEM_y,DEM_z); colormap(gca,jet); colorbar;
    set(gca,'ydir','normal','clim',[0 50]); hold on;
    set(gca,'xlim',[min(SL(berg_no).final.x)-150 max(SL(berg_no).final.x)+150],'ylim',[min(SL(berg_no).final.y)-150 max(SL(berg_no).final.y)+150]);
    plot(SL(berg_no).final.x,SL(berg_no).final.y,'-k','linewidth',2); hold on; plot(SL(berg_no).final.x,SL(berg_no).final.y,'--w','linewidth',2); hold on;
    plot(SL(berg_no).final.x_perim,SL(berg_no).final.y_perim,'-w','linewidth',2); hold on; plot(SL(berg_no).final.x_perim,SL(berg_no).final.y_perim,'--m','linewidth',2); hold on;
    
    %use the mask to extract size & shape info
    disp('Pulling plan-view info');
    for k = 1:size(x,1)-1
        dist(k) = ((x(k)-x(k+1))^2+(y(k)-y(k+1))^2)^(1/2);
    end
    dist(end+1) = ((x(end)-x(1))^2+(y(end)-y(1))^2)^(1/2);
    shapes_f = regionprops(iceberg_IMmask,'Area','MinorAxisLength','MajorAxisLength');
    lmin_f = abs(B.x(2)-B.x(1))*shapes_f.MinorAxisLength;
    lmax_f = abs(B.x(2)-B.x(1))*shapes_f.MajorAxisLength;
    SL(berg_no).final.width_min = lmin_f; SL(berg_no).final.width_max = lmax_f;
    base_area = shapes_f.Area.*im2_pixel_area; SL(berg_no).final.SA = base_area;
    SL(berg_no).final.z_map = single(DEM_z_masked); SL(berg_no).final.z_map(SL(berg_no).final.z_map == 0) = NaN; SL(berg_no).final.z_map(SL(berg_no).final.z_map<0) = NaN;
    zf = SL(berg_no).final.z_map(~isnan(SL(berg_no).final.z_map));
    SL(berg_no).final.z_median = nanmedian(zf); SL(berg_no).final.z_mad = mad(zf,1);
    perimeter = sum(dist);
    SL(berg_no).final.radius = perimeter/(2*pi);
    clear iceberg_IMmask;
    
    %iteratively estimate bulk density
    if SL(i).orientation == 0 %berg is a fragment or overturned so assume firn is totally saturated
        for k = 1:3 %k=1 is the best guess, k=[2,3] constraints uncertainties
            j=1;
            Hberg(j) = (rho_sw/(rho_sw-rho_i))*SL(i).final.z_median;
            Hberg_err(j) = abs(Hberg(j)).*sqrt((rho_sw_err/rho_sw)^2 + ((1.4826*SL(i).final.z_mad)/SL(i).final.z_median)^2 + ((rho_sw_err^2 + rho_i_err^2)/(rho_sw-rho_i)^2));
            wet_ref = find(density_z<=density.eightthir,1,'last');
            if Hberg(j) > density_z(wet_ref)
                rho_prof = [wetdensity_profile(1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [minwetdensity_profile(1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxwetdensity_profile(1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
            else
                rho_prof = [wetdensity_profile(1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [minwetdensity_profile(1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxwetdensity_profile(1:ceil(Hberg(j))+1)];
            end
            rho_f(j) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
            rho_f_err(j) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]));
            clear rho_prof*;

            while j
                Hberg(j+1) = (rho_sw/(rho_sw-rho_f(j)))*SL(i).final.z_median;
                Hberg_err(j+1) = abs(Hberg(j+1)).*sqrt(((abs(rho_sw/(rho_sw-rho_f(j)))*sqrt((rho_sw_err/rho_sw)^2 + (sqrt(rho_sw_err^2+rho_f_err(j)^2)/(rho_sw-rho_f(j)))^2))/(rho_sw/(rho_sw-rho_f(j))))^2 + ((1.4826*SL(i).final.z_mad)/SL(i).final.z_median)^2);
                if Hberg(j+1) > density_z(wet_ref)
                    rho_prof = [wetdensity_profile(1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [minwetdensity_profile(1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxwetdensity_profile(1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                else
                    rho_prof = [wetdensity_profile(1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [minwetdensity_profile(1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxwetdensity_profile(1:ceil(Hberg(j+1))+1)];
                end
                rho_f(j+1) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
                rho_f_err(j+1) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]));
                clear rho_prof*;

                if abs(rho_f(j+1)-rho_f(j)) < 0.25*rho_f_err(j+1)
                    if k == 1
                        SL(i).density_type = 'wet';
                        SL(i).final.density = rho_f(j+1);
                        %propagated density uncertainty
                        SL(i).final.density_uncert = rho_f_err(j);
                    else
                        %density uncertainty using FAC uncertainty bounds to iterate thickness & rho_f estimates
                        SL(i).final_range.density(k-1) = rho_f(j+1);
                    end
                    clear Hberg rho_f* dry_ref wet_ref;
                    break
                else
                    j = j+1;
                end
%                 clear Hice Vice;
            end
        end
    else %berg is upright so only wet the firn below the waterline
        for k = 1:3 %k=1 is the best guess, k=[2,3] constraints uncertainties
            j=1;
            Hberg(j) = (rho_sw/(rho_sw-rho_i))*SL(i).final.z_median;
            Hberg_err(j) = abs(Hberg(j)).*sqrt((rho_sw_err/rho_sw)^2 + ((1.4826*SL(i).final.z_mad)/SL(i).final.z_median)^2 + ((rho_sw_err^2 + rho_i_err^2)/(rho_sw-rho_i)^2));
            dry_ref = find(density_z<=SL(i).final.z_median,1,'last'); wet_ref = find(density_z<=density.eightthir,1,'last');
            if Hberg(j) > density_z(wet_ref)
                rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
            else
                rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:ceil(Hberg(j))+1)];
                rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:ceil(Hberg(j))+1)];
            end
            rho_f(j) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
            rho_f_err(j) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]));
            clear rho_prof*;

            while j
                Hberg(j+1) = (rho_sw/(rho_sw-rho_f(j)))*SL(i).final.z_median;
                Hberg_err(j+1) = abs(Hberg(j+1)).*sqrt(((abs(rho_sw/(rho_sw-rho_f(j)))*sqrt((rho_sw_err/rho_sw)^2 + (sqrt(rho_sw_err^2+rho_f_err(j)^2)/(rho_sw-rho_f(j)))^2))/(rho_sw/(rho_sw-rho_f(j))))^2 + ((1.4826*SL(i).final.z_mad)/SL(i).final.z_median)^2);
                if Hberg(j+1) > density_z(wet_ref)
                    rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
                else
                    rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(1,:) = [mindensity_profile(1:dry_ref) minwetdensity_profile(dry_ref+1:ceil(Hberg(j+1))+1)];
                    rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) maxwetdensity_profile(dry_ref+1:ceil(Hberg(j+1))+1)];
                end
                rho_f(j+1) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
                rho_f_err(j+1) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]));
                clear rho_prof*;

                if abs(rho_f(j+1)-rho_f(j)) < 0.25*rho_f_err(j+1)
                    if k == 1
                        SL(i).density_type = 'dry';
                        SL(i).final.density = rho_f(j+1);
                        %propagated density uncertainty
                        SL(i).final.density_uncert = rho_f_err(j);
                    else
                        %density uncertainty using FAC uncertainty bounds to iterate thickness & rho_f estimates
                        SL(i).final_range.density(k-1) = rho_f(j+1);
                    end
                    clear Hberg rho_f* dry_ref wet_ref;
                    break
                else
                    j = j+1;
                end
%                 clear Hice Vice;
            end
        end
    end
    
    %estimate the iceberg depth & submerged area using densities from the Ligtenberg FDM
    disp('Estimating thickness & volume');
    if SL(berg_no).orientation == 0 %overturned/fragment
        rho_f = sort([rho_i SL(berg_no).final.density]); %figure out if bubble-free ice or water-saturated firn has a lower density
        for l = 1:2
            draft(l) = (rho_f(l)/(rho_sw-rho_f(l)))*SL(berg_no).final.z_median;
            Hberg(l) = (rho_sw/(rho_sw-rho_f(l)))*SL(berg_no).final.z_median;
        end
    else
        rho_f = sort(SL(berg_no).final_range.density); %assume density constrained by model FAC range
        for l = 1:2
            draft(l) = (rho_f(l)/(rho_sw-rho_f(l)))*SL(berg_no).final.z_median;
            Hberg(l) = (rho_sw/(rho_sw-rho_f(l)))*SL(berg_no).final.z_median;
        end
    end
    SL(berg_no).final_range.draft = sort(draft);
    lat_area = sort(draft).*sum(dist); SL(berg_no).final_range.LA = lat_area;
    area = base_area*ones(size(lat_area)) + lat_area; SL(berg_no).final_range.TA = area;
    SL(berg_no).final.V = DEM2_pixel_area*(rho_sw/(rho_sw-SL(berg_no).final.density))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(berg_no).final.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(berg_no).final_range.V(1) = DEM2_pixel_area*(rho_sw/(rho_sw-min(rho_f)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(berg_no).final.z_median*sum(sum(isnan(DEM_z_masked))));
    SL(berg_no).final_range.V(2) = DEM2_pixel_area*(rho_sw/(rho_sw-max(rho_f)))*(sum(DEM_z_masked(~isnan(DEM_z_masked)))+SL(berg_no).final.z_median*sum(sum(isnan(DEM_z_masked))));
    
    %save the iceberg outline
    S.Geometry = 'Polygon';
    S.BoundingBox = [min(x) min(y); max(x) max(y)];
    S.X = double(x'); S.Y = double(y');
    S.Name = ['iceberg',num2str(berg_number)];
    shapefile_name = ['WV_',num2str(tf),'_icebergshape',num2str(berg_number)];
    cd([dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/']);
    shapewrite(S,shapefile_name);
    copyfile([dir_code,'general/antarctic_PSprojection.prj'],[dir_output,DEM1.time,'-',DEM2.time,'/iceberg_shapes/',shapefile_name,'.prj']);
    cd([dir_output,DEM1.time,'-',DEM2.time,'/']);
    clear S;
    
    %calculate mean surface area for the exposed face (SA) &
    %submerged (TA) portions of the iceberg & associated temporal uncertainty
    SL(berg_no).mean.z = (SL(berg_no).initial.z_median+SL(berg_no).final.z_median)/2;
    SL(berg_no).mean.SA = (SL(berg_no).initial.SA+SL(berg_no).final.SA)/2;
    SL(berg_no).change.SA = abs(SL(berg_no).initial.SA-SL(berg_no).final.SA)/2;
    SL(berg_no).mean.TA = (nanmean(SL(berg_no).initial_range.TA)+nanmean(SL(berg_no).final_range.TA))/2;
    SL(berg_no).uncert.TA = max([(SL(berg_no).initial_range.TA(2)-SL(berg_no).initial_range.TA(1))/2 (SL(berg_no).final_range.TA(2)-SL(berg_no).final_range.TA(1))/2]);
    SL(berg_no).change.TA = abs(nanmean(SL(berg_no).initial_range.TA)-nanmean(SL(berg_no).final_range.TA))/2;
    SL(berg_no).mean.V = (nanmean(SL(berg_no).initial_range.V)+nanmean(SL(berg_no).final_range.V))/2;
    SL(berg_no).uncert.V = max([(SL(berg_no).initial_range.V(2)-SL(berg_no).initial_range.V(1))/2 (SL(berg_no).final_range.V(2)-SL(berg_no).final_range.V(1))/2]);
    SL(berg_no).change.V = abs(nanmean(SL(berg_no).initial_range.V)-nanmean(SL(berg_no).final_range.V))/2;
    
    %save the data
    disp('Saving data');
    save([dir_iceberg,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat'],'SL','-v7.3');
    clear DEM_z_masked z idx idy x y iceberg_IMmask iceberg_DEMmask rho_f draft Hberg lat_area base_area dist area_mask masked_pixels image_xgrid image_ygrid shapes lmin* lmax* perimeter;
    clear x1 x2 y1 y2 xlims ylims xmin xmax ymin ymax DEM_x DEM_y DEM_z;
    close all; drawnow;
    
    
    %cd to the melt rate file & load each file sequentially
    cd(dir_iceberg);
    if i < 10
        load_file = ['load iceberg',num2str(0),num2str(i),'_dz.mat']; eval(load_file);
    else
        load_file = ['load iceberg',num2str(i),'_dz.mat']; eval(load_file);
    end
    
    %estimate the density
    if SL(berg_no).orientation == 0
        rho_f = nanmean([SL(berg_no).initial.density SL(berg_no).final.density]);
        rho_f_range = [rho_i SL(berg_no).initial_range.density SL(berg_no).final_range.density];
        rho_f_err = [min(rho_f_range-rho_f) max(rho_f_range-rho_f)];
    else
        rho_f = nanmean([SL(berg_no).initial.density SL(berg_no).final.density]);
        rho_f_range = [SL(berg_no).initial_range.density SL(berg_no).final_range.density];
        rho_f_err = [min(rho_f_range-rho_f) max(rho_f_range-rho_f)];
    end
    
    %calculate the ice thickness (change to H = 3*(V/SA) if assuming the submerged shape is a cone!)
    dt=sum(SL(berg_no).days);
    SL(berg_no).initial_range.H = (repmat(rho_sw,size(SL(berg_no).initial_range.density))./(repmat(rho_sw,size(SL(berg_no).initial_range.density))-SL(berg_no).initial_range.density)).*SL(berg_no).initial.z_median;
    SL(berg_no).final_range.H = (repmat(rho_sw,size(SL(berg_no).final_range.density))./(repmat(rho_sw,size(SL(berg_no).final_range.density))-SL(berg_no).final_range.density)).*SL(berg_no).final.z_median;
    SL(berg_no).mean.H = nanmean([SL(berg_no).initial_range.H SL(berg_no).final_range.H]);
    SL(berg_no).change.H = nanmean(SL(berg_no).initial_range.H) - nanmean(SL(berg_no).final_range.H);
    SL(berg_no).mean.draft = nanmean([SL(berg_no).initial_range.draft SL(berg_no).final_range.draft]);
    SL(berg_no).change.draft = nanmean(SL(berg_no).initial_range.draft) - nanmean(SL(berg_no).final_range.draft);
    
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
    SL(berg_no).mean.dz = nanmean(vals); SL(berg_no).uncert.dz = nanstd(vals);
    dZ_mean = (rho_sw/(rho_sw-nanmean([SL(berg_no).initial.density SL(berg_no).final.density])))*SL(berg_no).mean.dz;
    dZ_stdev = (rho_sw/(rho_sw-nanmean([SL(berg_no).initial.density SL(berg_no).final.density])))*SL(berg_no).uncert.dz;
    
    %correct for SMB loss
    dH_SMBadjust_mean = dZ_mean + SL(berg_no).SMB; %SL(berg_no).SMB=-surfmelt
    
    %correct for creep thinning
    rf = SL(berg_no).ratefactor^(-1/3); %Pa s^1/3
    creep = ((-1/(2*sqrt(3)))*((rho_f*9.81*SL(berg_no).mean.z)/(2*sqrt(3)))^3*(1-(rho_f/rho_sw))^3)/(rf^3); %creep thinning rate (1/s)
    SL(berg_no).creep_dz = (SL(berg_no).mean.H*creep*(dt*31536000));
    dH_submelt = dH_SMBadjust_mean + SL(berg_no).creep_dz; %integrate creep over the ice thickness & over the time period
    
    %dH uncertainty sources
    %S1) Bulk density (use range in density estimates, including ice density if overturned)
    %S2) Ocean water density (use 1026 +/- 2 kg/m^3)
    %S3) Runoff uncertainty (from RACMO output)
    %R1) DEM uncertainty (random errors estimated as 2.9m)
    %R2) Translocation & Rotation (User) error (quantify by repeating procedure 10 times for each glacier)
    
    %%quantify potential bias from systematic errors
    %uncertainty from S1 & S2
    rho_added_err = sqrt(rho_f_err.^2 + repmat(rho_sw_err,size(rho_f_err)).^2);
    rho_conversion = rho_sw/(rho_sw-rho_f);
    rho_conversion_err = abs(rho_conversion)*sqrt(repmat((rho_sw_err/rho_sw),size(rho_f_err)).^2 + (rho_added_err./repmat(rho_sw-rho_f,size(rho_f_err))).^2);
    
    %uncertainty from S3
    dz_SMB_err = 0.3*abs(SL(berg_no).SMB);
    
    %uncertainty from R1
    pixelno_o = SL(berg_no).initial.SA./DEM1_pixel_area; pixelno_f = SL(berg_no).final.SA./DEM2_pixel_area;
    dz_randerr = sqrt((2.9^2+2.9^2)./nanmean([pixelno_o pixelno_f])); %(2.9^2+2.9^2) is to account for DEM differencing over time
    
    %uncertainty from R2
    for k = 1:10
        dzs(k) = IB(k).dz.local_adjust.mean;
    end
    dz_user_err = nanstd(dzs);
    
    %total, systematic, and random errors
    SL(berg_no).uncert.h_SMB = dz_SMB_err;
    SL(berg_no).uncert.h_rand = dz_randerr;
    SL(berg_no).uncert.h_user = dz_user_err;
    SL(berg_no).uncert.hH_conversion = rho_conversion_err;
    dH_total_err = abs(dH_submelt).*sqrt(repmat(((dz_SMB_err^2 + dz_randerr^2 + dz_user_err^2)/SL(berg_no).mean.dz^2),size(rho_conversion_err)) + (rho_conversion_err./repmat(rho_conversion,size(rho_conversion_err))).^2);
    
    %convert all 'dHs', which are change in thickness as estimated using the elevation change,
    %(& errors) to changes in volume so that the iceberg melt rate
    %(dH/dt) can be calculated by dividing the volume change by the total
    %submerged iceberg area
    dV_mean = SL(berg_no).mean.SA*dH_submelt;
    dV_stdev = abs(dV_mean).*sqrt((dZ_stdev./dH_submelt)^2);
    dV_total_err = abs(dV_mean).*sqrt((dH_total_err./repmat(dH_submelt,size(dH_total_err))).^2 + repmat((SL(berg_no).change.SA/2)./SL(berg_no).mean.SA,size(dH_total_err)).^2);
    
    %convert to a melt rate (volume change with time & thickness change
    %over the entire submerged face)
    dVdt_mean = dV_mean/dt;  dVdt_stdev = dV_stdev/dt;
    dVdt_total_err = dV_total_err./dt;
    dHdt_mean = dVdt_mean/SL(berg_no).mean.TA;
    dHdt_stdev = abs(dHdt_mean).*sqrt((dVdt_stdev./dVdt_mean)^2+(SL(berg_no).uncert.TA./SL(berg_no).mean.TA)^2);
    dHdt_total_err = abs(dHdt_mean).*sqrt((dVdt_total_err./repmat(dVdt_mean,size(dVdt_total_err))).^2 + repmat(SL(berg_no).uncert.TA./SL(berg_no).mean.TA,size(dVdt_total_err)).^2);
    
    %save volume fluxes and melt rates to the structure
    SL(berg_no).mean.dVdt = dVdt_mean;
    SL(berg_no).range.dVdt(1) = dVdt_mean-3*dVdt_stdev; SL(berg_no).range.dVdt(2) = dVdt_mean+3*dVdt_stdev;
    SL(berg_no).uncert.dVdt = dVdt_total_err;
    SL(berg_no).mean.dHdt = dHdt_mean;
    SL(berg_no).range.dHdt(1) = dHdt_mean-3*dHdt_stdev; SL(berg_no).range.dHdt(2) = dHdt_mean+3*dHdt_stdev;
    SL(berg_no).uncert.dHdt = dHdt_total_err;
    
    %if overlapping bedrock regions were present, add the uniform offset that
    %would be applied from bedrock & tidal differences to the structure
%     if ~isnan(IB(1).dz.br_tide_adjust.mean)
%         SL(berg_no).elpsd_adjust_o = IB(1).elpsd_adjust_o; SL(berg_no).elpsd_bt_adjust_f = IB(1).elpsd_adjust_f+IB(1).dbedrock_f+IB(1).dtide_f;
%   end
    
end

%save the compiled data
cd(dir_iceberg);
disp('Saving melt rates');
save([dir_output,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat'],'SL','-v7.3');
disp('Melt rates resaved!');
close all; drawnow;
clear A B Y Z;

%plot to figure-out which icebergs need to be removed
dVdt = []; Asub = []; H = []; m = []; coreg_zo = []; coreg_zf = []; berg_ref =[];
for i = 1:length(SL)
    if SL(i).mean.dVdt > 0 && ~isempty(SL(i).mean.TA)
        dVdt = [dVdt SL(i).mean.dVdt];
        Asub = [Asub SL(i).mean.TA];
        H = [H SL(i).mean.H];
        m = [m SL(i).mean.dHdt];
        coreg_zo = [coreg_zo SL(i).initial.coreg_z]; coreg_zf = [coreg_zf SL(i).final.coreg_z];
        berg_ref = [berg_ref; SL(i).name(end-1:end)];
    end
end
figure; set(gcf,'position',[100 500 1500 600]);
subplot(1,3,1);
plot(Asub,dVdt,'ok','markersize',24,'markerfacecolor','w'); hold on;
set(gca,'fontsize',20); xlabel('Submerged area (m^2)','fontsize',20); ylabel('Meltwater flux (m^3/d)','fontsize',20);
for i = 1:length(m)
    text(double(Asub(i))-0.03e5,double(dVdt(i)),berg_ref(i,:))
end
grid on;
subplot(1,3,2);
plot(H,m,'ok','markersize',24,'markerfacecolor','w'); hold on;
set(gca,'fontsize',20); xlabel('Average iceberg thickness (m)','fontsize',20); ylabel('Melt rate (m/d)','fontsize',20);
for i = 1:length(m)
    text(double(H(i))-2,double(m(i)),berg_ref(i,:))
end
grid on;
subplot(1,3,3);
plot(coreg_zo-coreg_zf,dVdt,'ok','markersize',24,'markerfacecolor','w'); hold on;
set(gca,'fontsize',20); xlabel('\Delta sea-level adjustment (m)','fontsize',20); ylabel('Meltwater flux (m^3/d)','fontsize',20);
for i = 1:length(m)
    text(double(coreg_zo(i)-coreg_zf(i)),double(dVdt(i)),berg_ref(i,:))
end
grid on;
disp('Iceberg meltwater flux should increase linearly with submerged area');
disp('Iceberg melt rates should increase with thickness');

%replace bad melt rate and area estimates with empty brackets
prompt = 'Do you want to remove any icebergs from the analysis  (y/n)?';
str = input(prompt,'s');
if strmatch(str,'y')==1
    disp(['Already identified bad icebergs as #s ',num2str(bad_bergs)]);
    disp('Specify the icebergs to remove as "bad_bergs=[A B C etc]; return"');
    keyboard
    for i =bad_bergs
        SL(i).mean.TA = [];
        SL(i).mean.dHdt = [];
    end
end
save_file = ['save(''',region_abbrev,'_iceberg_melt.mat''',',','''SL''',',','''-v7.3''',')']; eval(save_file);

end