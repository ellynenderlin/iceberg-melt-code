%extract information on firn air content for the region
clear all; close all;
region_name = 'LarsenB'; region_abbrev = 'LB';
root_dir = '/Users/enderlin/Antarctic_icebergs/'; 


%% extract terminus elevation data & delineate areas occupied by icebergs
cd_root_dir = ['cd ',root_dir]; eval(cd_root_dir);
disp('loading the low-res RAMP DEM to fill data gaps');
[A,R] = geotiffread('Antarctic_RAMP_DEM.tif');
RAMP.x = R.XWorldLimits(1)+0.5*R.CellExtentInWorldX:R.CellExtentInWorldX:R.XWorldLimits(2)-0.5*R.CellExtentInWorldX;
RAMP.y = R.YWorldLimits(2)-0.5*R.CellExtentInWorldY:-R.CellExtentInWorldY:R.YWorldLimits(1)+0.5*R.CellExtentInWorldY;
RAMP.z = single(A); RAMP.z(RAMP.z == -9999) = NaN;

%read-in all the elevation datasets
disp('reading-in all DEMs for the region & associated iceberg coordinates');
cd_region_dir = ['cd ',root_dir,'/',region_name]; eval(cd_region_dir);
DEMdirs = dir('20*');
for i = 1:length(DEMdirs)
    disp(['Looping through date directory #',num2str(i),' of ',num2str(length(DEMdirs))]);
    cd_DEM_dir = ['cd ',DEMdirs(i).name]; eval(cd_DEM_dir);
    DEMs = dir('*DEM.mat');
    for j = 1:length(DEMs)
        load_DEM = ['load ',DEMs(j).name]; eval(load_DEM);
        w=who;
        if ~isempty(strmatch('Z',w))
            D(i).x=Z.x; D(i).y = Z.y;
            D(i).z = Z.z_elpsd_adjust;
        else
            D(i).x=Y.x; D(i).y = Y.y;
            D(i).z = Y.z_elpsd_adjust;
        end
        
        load_image = ['load ',DEMs(j).name(1:end-7),'orthoimage.mat']; eval(load_image);
        I(j).x = IM.x; I(j).y = IM.y; I(j).z = IM.z;
    end
    [D(i).xgrid, D(i).ygrid] = meshgrid(D(i).x,D(i).y);
    
    %load iceberg coordinates
    cd iceberg_data
    if i == 1; iceberg_x = []; iceberg_y = []; end
    iceberg_coords = dir('*PScoords.txt');
    for j = 1:length(iceberg_coords)
        coords = dlmread(iceberg_coords(j).name);
        iceberg_x = [iceberg_x; nanmean([coords(2) coords(4)]);];
        iceberg_y = [iceberg_y; nanmean([coords(1) coords(3)]);];
    end
    
    cd ../..
end

%interpolate elevations to a standard grid
xlims=[]; ylims=[];
for i = 1:length(D)
    xlims(i,:) = [min(D(i).x) max(D(i).x)]; ylims(i,:) = [min(D(i).y) max(D(i).y)];
end
RAMPxref(1) = find(RAMP.x<=min(min(xlims)),1,'last'); RAMPxref(2) = find(RAMP.x>=max(max(xlims)),1,'first');
RAMPyref(2) = find(RAMP.y>=min(min(ylims)),1,'last'); RAMPyref(1) = find(RAMP.y<=max(max(ylims)),1,'first');
ramp.x = RAMP.x(min(RAMPxref):max(RAMPxref)); ramp.y = RAMP.y(min(RAMPyref):max(RAMPyref)); 
ramp.z = RAMP.z(min(RAMPyref):max(RAMPyref),min(RAMPxref):max(RAMPxref));
h.x = min(min(xlims)):5:max(max(xlims)); h.y = min(min(ylims)):5:max(max(ylims));
[h.xgrid,h.ygrid] = meshgrid(h.x,h.y);
for i = 1:length(D)
   D(i).z_interp = interp2(D(i).xgrid,D(i).ygrid,D(i).z,h.xgrid,h.ygrid); 
   z(:,:,i) = D(i).z_interp;
end
h.z = nanmean(z,3);

%fill-in NaN regions with low-res RAMP DEM to provide regional context
[rampxgrid,rampygrid] = meshgrid(ramp.x,ramp.y);
ramp.z_interp = interp2(rampxgrid,rampygrid,ramp.z,h.xgrid,h.ygrid); 
ramp.z_interp = ramp.z_interp - nanmean(min(ramp.z_interp));
h.z_filled = h.z;
h.z_filled(isnan(h.z)) = ramp.z_interp(isnan(h.z));

%plot average DEM & map glacier termini and the extent of calved icebergs (if multiple termini are present)
figure; set(gcf,'position',[50 50 800 800]);
cmap = colormap(jet(10001));
imagesc(h.x,h.y,h.z_filled); axis xy equal; colormap(gca,cmap); hold on;
set(gca,'clim',[floor(nanmean(min(h.z)))-1 floor(nanmean(min(h.z)))+99]);
drawnow;
%plot the coordinates of the icebergs you've identified for all dates
for i = 1:length(iceberg_x)
   plot(iceberg_x(i),iceberg_y(i),'pk','markerfacecolor','y'); hold on; 
end
drawnow;
%plot a border around the high-res DEM
for i = 1:length(h.y)
    if sum(~isnan(h.z(i,:))) >=1
    plot(h.x(find(~isnan(h.z(i,:)),1,'first')),h.y(i),'.w'); hold on;
    plot(h.x(find(~isnan(h.z(i,:)),1,'last')),h.y(i),'.w'); hold on;
    end
end
drawnow;
disp('Extract terminus elevations & delineate the area containing icebergs calved from each glacier');
j=1;
while j
    disp('map the glacier terminus');
    disp('click the UL & LR corners of a box surrounding the glacier terminus & calved icebergs to zoom');
    a = ginput(2);
    set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
    disp('draw a polygon around the area within ~1km from the terminus');
    [tz,tx,ty] = roipoly; tmasked = tz.*h.z; tmasked(tz==0) = NaN;
    term_mask(j).x = tx; term_mask(j).y = ty; term_mask(j).mask = tz;
    term_mask(j).zmasked = tmasked;
    term_mask(j).zmean = nanmean(term_mask(j).zmasked(~isnan(term_mask(j).zmasked)));
    term_mask(j).zstd = nanstd(term_mask(j).zmasked(~isnan(term_mask(j).zmasked)));
    term_mask(j).zmedian = nanmedian(term_mask(j).zmasked(~isnan(term_mask(j).zmasked)));
    term_mask(j).zmad = mad(term_mask(j).zmasked(~isnan(term_mask(j).zmasked)),1);
    disp('draw a polygon to roughly delineate the area containing icebergs calved from the glacier');
    [~,ix,iy] = roipoly; iceberg_mask(j).x = ix; iceberg_mask(j).y = iy;
    plot(iceberg_mask(j).x,iceberg_mask(j).y,'-m'); hold on;
    set(gca,'xlim',[min(h.x) max(h.x)],'ylim',[min(h.y) max(h.y)]);
    clear tx ty tz ix iy;
    
    prompt = 'Is there another source glacier from which icebergs calved (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        disp('complete the same steps for the other glacier');
        j=j+1; 
    else
        break
    end
    
end
save_data = ['save ',region_name,'_calving_data.mat h term_mask iceberg_mask']; eval(save_data);
close all; drawnow;
disp('Saved terminus elevation & iceberg extent info');


%% determine whether each glacier terminus is floating or grounded
cd_region_dir = ['cd ',root_dir,'/',region_name]; eval(cd_region_dir);
load_data = ['load ',region_name,'_calving_data.mat']; eval(load_data);
fig1 = figure; set(gcf,'position',[50 50 800 800]);
cmap = colormap(jet(10001));
imagesc(h.x,h.y,h.z_filled); axis xy equal; colormap(gca,cmap); hold on;
set(gca,'clim',[floor(nanmean(min(h.z)))-1 floor(nanmean(min(h.z)))+99]);
drawnow;

%read in the Depoorter et al. (2013) Antarctic grounding line dataset
cd_root_dir = ['cd ',root_dir]; eval(cd_root_dir);
GL=dlmread('Antarctic_GL-coords.txt');
for i = 1:length(GL)
    [GLx(i),GLy(i)] = wgs2ps(GL(i,2),GL(i,1),'StandardParallel',-71,'StandardMeridian',0);
end

%determine if the glacier is floating
fig2 = figure; set(gcf,'position',[800 50+600-(200*length(I)/2) 700 300*length(I)/2]);
graymap = colormap(gray(10001));
for k = 1:length(I)
    subp(k) = subplot(ceil(length(I)/2),2,k);
    imagesc(I(k).x,I(k).y,I(k).z); axis xy equal; colormap(gca,graymap);
end
drawnow;
for j = 1:length(term_mask)
    
    %locate the Depoorter grounding line
    dx = GLx-nanmean(term_mask(j).x); dy = GLy-nanmean(term_mask(j).y);
    diff = sqrt(dx.^2 + dy.^2);
    GLdist(j) = min(diff);
    GLref(j) = find(diff==min(diff));
    clear dx dy diff;
    
    %plot a portion of the grounding line on the figure
    figure(fig1);
    plot(GLx(GLref(j)-10:GLref(j)+10),GLy(GLref(j)-10:GLref(j)+10),'.k','markersize',12); hold on
    drawnow;
    
    %extract a near-terminus along-flow elevation profile if possible
    prompt = 'Is there high-res elevation data over the terminus (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        disp('trace a flowline from the terminus up the glacier');
        [px,py,pz] = improfile;
        proffig = figure; set(proffig,'position',[700 50 600 300]);
        plot(sqrt((px-px(1)).^2 + (py-py(1)).^2),smooth(pz,21),'-k'); hold on;
        clear px py pz
    end
    drawnow;
    
    %zoom in on the icebergs 
    figure(fig2);
    for k = 1:length(I)
    set(subp(k),'xlim',[min(iceberg_mask(j).x)-1000 max(iceberg_mask(j).x)+1000],...
        'ylim',[min(iceberg_mask(j).y)-1000 max(iceberg_mask(j).y)+1000]);
    end
    drawnow;
    
    %based on the proximity of the Depoorter grounding line, the glacier surface slope, 
    %& whether you saw any upright icebergs in the images, decide if you think the glacier is grounded
    disp('General characteristics of floating termini (in order of importance/trustworthiness):');
    disp('1) Upright icebergs (zoom in on subplot images to look at icebergs)');
    disp('2) Nearly flat surface slope near the terminus');
    disp('3) Close proximity to Depoorter grounding line');
    prompt = 'Is the terminus floating (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        term_mask(j).floating = 1;
    else
        term_mask(j).floating = 0;
    end

end
cd_region_dir = ['cd ',root_dir,'/',region_name]; eval(cd_region_dir);
save_data = ['save ',region_name,'_calving_data.mat h term_mask iceberg_mask']; eval(save_data);
close all; drawnow;


%% read in firn density data & come up with uncertainty estimates
cd_to_root_dir = ['cd ',root_dir]; eval(cd_to_root_dir); cd FDM_Antarctica
%FAC = firn air content (H_observed-FAC = H_i)
prompt = 'Are you looking at icebergs along the Antarctic Peninsula (y/n)?';
str = input(prompt,'s');
if strmatch(str,'y')==1
    firn_lat = ncread('FDM_FirnAir_XPEN055_1979-2016.nc','lat');
    firn_lon = ncread('FDM_FirnAir_XPEN055_1979-2016.nc','lon');
    FAC = ncread('FDM_FirnAir_XPEN055_1979-2016.nc','FirnAir'); 
    firnu_lat = ncread('FDM_FirnAir-uncertainty_XPEN055_1979-2016.nc','lat'); 
    firnu_lon = ncread('FDM_FirnAir-uncertainty_XPEN055_1979-2016.nc','lon'); 
    FACu = ncread('FDM_FirnAir-uncertainty_XPEN055_1979-2016.nc','Total_unc'); 
    firn_time = ncread('FDM_FirnAir_XPEN055_1979-2016.nc','time');
else
    firn_lat = ncread('FDM_FirnAir_ANT27_1979-2016.nc','lat');
    firn_lon = ncread('FDM_FirnAir_ANT27_1979-2016.nc','lon');
    FAC = ncread('FDM_FirnAir_ANT27_1979-2016.nc','FirnAir'); 
    firnu_lat = ncread('FDM_FirnAir-uncertainty_ANT27_1979-2016.nc','lat'); 
    firnu_lon = ncread('FDM_FirnAir-uncertainty_ANT27_1979-2016.nc','lon'); 
    FACu = ncread('FDM_FirnAir-uncertainty_ANT27_1979-2016.nc','Total_unc'); 
    firn_time = ncread('FDM_FirnAir_ANT27_1979-2016.nc','time');
end
for i = 1:size(firn_lat,1)
    for j = 1:size(firn_lat,2)
        [firn_x(i,j),firn_y(i,j)] = wgs2ps(firn_lon(i,j),firn_lat(i,j),'StandardParallel',-71,'StandardMeridian',0);
    end
end
clear *diff*;
lat_diff = abs(nanmean(iceberg_y)*ones(size(firn_y)) - firn_y);
lon_diff = abs(nanmean(iceberg_x)*ones(size(firn_x)) - firn_x);
diff_map = sqrt(lat_diff.^2+lon_diff.^2);
firn_ref = find(diff_map==min(min(diff_map)));
[firny firnx] = ind2sub(size(nanmean(FAC(:,:,size(FAC,3)-37:size(FAC,3)),3)),firn_ref);
disp(['FAC x-reference = ',num2str(firnx),' & y-reference = ',num2str(firny)]);

%adjust FAC reference grid cell if necessary
disp('Adjust coordinates (if necessary) to extract surface melt estimates');
figure; set(gcf,'position',[500 100 900 900]);
imagesc(nanmean(FAC,3)); colormap(jet(10001)); hold on; set(gca,'ydir','reverse'); 
cbar = colorbar; set(get(cbar,'ylabel'),'string','FAC (m) ');
plot(firnx,firny,'ok','markerfacecolor','k','markeredgecolor','w','markersize',20); hold on;
plot(firnx,firny,'xw','markersize',20); hold on;
set(gca,'xlim',[firnx-10 firnx+10],'ylim',[firny-10 firny+10]);
prompt = 'Modify coordinates if the marker is in a region with no data. Do the coordinates need to be modified (y/n)?';
str = input(prompt,'s');
if strmatch(str,'y')==1
    disp('To change coordinates, identify new coordinates using tick marks and type firnx=XX; firny=YY; dbcont ');
    keyboard
end
close all; drawnow;

%interpolate the FAC uncertainty map to the FAC grid
for i = 1:size(firnu_lat,1)
    for j = 1:size(firnu_lat,2)
        [firnu_x(i,j),firnu_y(i,j)] = wgs2ps(firnu_lon(i,j),firnu_lat(i,j),'StandardParallel',-71,'StandardMeridian',0);
    end
end
facu_x=reshape(firnu_x,[(size(FACu,1)*size(FACu,2)),1]);
facu_y=reshape(firnu_y,[(size(FACu,1)*size(FACu,2)),1]);
facu=reshape(FACu,[(size(FACu,1)*size(FACu,2)),1]);
F=scatteredInterpolant(facu_x,facu_y,facu);
FACu_interp = F(firn_x,firn_y);

%pull the annual avg FAC & extrapolate density depth profiles to estimate
%uncertainties in iceberg density
firnair.xref=firnx; firnair.yref=firny;
firnair.mean = nanmean(FAC(firny,firnx,:)); firnair.std = nanstd(FAC(firny,firnx,:));
firnair.median = nanmedian(FAC(firny,firnx,:)); firnair.mad = mad(FAC(firny,firnx,:),1);
firnair.uncert = FACu_interp(firny,firnx);
density_lat = ncread('FDM_DepthDensLevels_ANT27_map.nc','lat');
density_lon = ncread('FDM_DepthDensLevels_ANT27_map.nc','lon');
sevenhundred = ncread('FDM_DepthDensLevels_ANT27_map.nc','z700');
sevenfifty = ncread('FDM_DepthDensLevels_ANT27_map.nc','z750');
eighthundred = ncread('FDM_DepthDensLevels_ANT27_map.nc','z800');
eightthirty = ncread('FDM_DepthDensLevels_ANT27_map.nc','z830');
for i = 1:size(density_lat,1)
    for j = 1:size(density_lat,2)
        [density_x(i,j),density_y(i,j)] = wgs2ps(density_lon(i,j),density_lat(i,j),'StandardParallel',-71,'StandardMeridian',0);
    end
end
clear *diff*;
lat_diff = abs(nanmean(iceberg_y)*ones(size(density_y)) - density_y);
lon_diff = abs(nanmean(iceberg_x)*ones(size(density_x)) - density_x);
diff_map = sqrt(lat_diff.^2+lon_diff.^2);
density_ref = find(diff_map==min(min(diff_map)));
[densityy,densityx] = ind2sub(size(eightthirty),density_ref);
disp(['Density profile x-reference = ',num2str(densityx),' & y-reference = ',num2str(densityy)]);
density.xref = densityx; density.yref = densityy;
density.sevhun = sevenhundred(density.yref,density.xref);
density.sevfif = sevenfifty(density.yref,density.xref);
density.eighthun = eighthundred(density.yref,density.xref);
density.eightthir = eightthirty(density.yref,density.xref);
density_levels = [700 750 800 830]; density_depths = [density.sevhun density.sevfif density.eighthun density.eightthir];
%fit a curve to the density-depth profile to estimate the depth of the base of the firn column
% [f,gof] = fit(density_levels',density_depths','smoothingspline');
% ex_density = fnxtr(f.p); %extrapolate outside the data domain using a second order polynomial
% density.nineseventeen = fnval(ex_density,917); %evaluate the extrapolation function at the bubble-free ice density
ft = fittype('917-(917-a)*exp(-x/b)');
[f,gof] = fit(density_depths',density_levels',ft,'StartPoint',[350 firnair.median]); %empirical, exponential density-depth relation from Schytt (1958) on Cuffey & Paterson p. 19
density.nineseventeen = -f.b*log(-(916.9-917)/(917-f.a)); %find depth where rho=916.9 (goes to infinity at 917)
density_levels = [density_levels 917]; density_depths = [density_depths density.nineseventeen];

%save the FAC & density data
cd_region_dir = ['cd ',root_dir,'/',region_name]; eval(cd_region_dir);
save_data = ['save ',region_name,'_calving_data.mat h term_mask iceberg_mask firnair density']; eval(save_data);
close all; drawnow;
