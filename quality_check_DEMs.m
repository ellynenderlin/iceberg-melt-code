%% Check DEM coverage & quality: delete if bad!
%Loops through files in a specified directory, identifies the DEMs, and
%loads an image of each DEM before prompting you to decide if you want to
%keep it. Bad DEMs and supporting files are deleted. Good DEMs will be
%reprojected so that elevations are orthometric (i.e., with respect to the
%geoid, with 0 m at sea level) instead of the native elllipsoidal
%elevations.

%Requirements:
% 1) Unique directory for each region.
% 2) Shapefile of the region outline in the directory
% 4) cmocean colormapping package
% 3) (optional but super helpful) Screenshot of image overlain with outline

%% initialize
clearvars; close all;

%specify site-specific parameters
glacier_abbrev = 'TKS'; %site name
root_dir = ['/Volumes/Growler/Greenland/melange/',glacier_abbrev,'/']; %update to your path 
% ROI = 'BB_APG.shp'; %update to the appropriate shapefile path & name (only need the name if in the directory with the DEMs)
ROI = ['/Users/Shared/Greenland/Bounding-Box-Shapefile/',glacier_abbrev,'/BB_',glacier_abbrev,'.shp']; %full path example

%pre-defined colormap
addpath('/Volumes/Growler/general-code/cmocean/','/Users/Shared/general-code/'); %UPDATE THIS PATH AS NEEDED!
elev_cmap = cmocean('thermal',1001); elev_cmap(1,:) = [1 1 1];

%% loop through files, trashing bad ones

%find files
cd(root_dir);
tifs = dir('SETSM*dem.tif'); %melange DEMs provided by PGC
S = shaperead(ROI); %fjord shapefile used by PGC to crop DEMs


%pull dates
for p = 1:length(tifs); DEMtif_dates(p,:) = tifs(p).name(19:26);end

%loop through & check files
for p = 1:length(tifs) %update the first number in here if you want to restart not at 1
    disp(['reformating DEM tif #',num2str(p),' of ',num2str(length(tifs))]);
    disp(DEMtif_dates(p,:));
    DEM_name = [glacier_abbrev,'-',DEMtif_dates(p,:),'_fjord-DEM.mat'];
    
    %load the geotiffs
    if contains(version,'R2019')
        [DEM,RDEM] = geotiffread(tifs(p).name); [IM,RIM] = geotiffread([tifs(p).name(1:end-7),'ortho.tif']); 
    else
        [DEM,RDEM] = readgeoraster(tifs(p).name); [IM,RIM] = readgeoraster([tifs(p).name(1:end-7),'ortho.tif']);
    end
    %create coordinate vectors
    if strcmp(RDEM.RasterInterpretation,'cells')
        Z.x = RDEM.XWorldLimits(1)+0.5*RDEM.CellExtentInWorldX:RDEM.CellExtentInWorldX:RDEM.XWorldLimits(2)-0.5*RDEM.CellExtentInWorldX;
        Z.y = RDEM.YWorldLimits(2)-0.5*RDEM.CellExtentInWorldY:-RDEM.CellExtentInWorldY:RDEM.YWorldLimits(1)+0.5*RDEM.CellExtentInWorldY;
        Y.x = RIM.XWorldLimits(1)+0.5*RIM.CellExtentInWorldX:RIM.CellExtentInWorldX:RIM.XWorldLimits(2)-0.5*RIM.CellExtentInWorldX;
        Y.y = RIM.YWorldLimits(2)-0.5*RIM.CellExtentInWorldY:-RIM.CellExtentInWorldY:RIM.YWorldLimits(1)+0.5*RIM.CellExtentInWorldY;
    else
        Z.x = RDEM.XWorldLimits(1):RDEM.CellExtentInWorldX:RDEM.XWorldLimits(2);
        Z.y = RDEM.YWorldLimits(2):-RDEM.CellExtentInWorldY:RDEM.YWorldLimits(1);
        Y.x = RIM.XWorldLimits(1):RIM.CellExtentInWorldX:RIM.XWorldLimits(2);
        Y.y = RIM.YWorldLimits(2):-RIM.CellExtentInWorldY:RIM.YWorldLimits(1);
    end
    Z.z.raw = double(DEM); Z.z.raw(Z.z.raw==-9999) = NaN;
    Y.z = double(IM); Y.z(Y.z == 32767) = NaN;
    zmap = Z.z.raw; zmap(Z.z.raw>200) = NaN; zmap(isinf(zmap)) = NaN;
    clear DEM RDEM IM RIM;
    
    %if all elevations are garbage, delete the file
    if sum(sum(~isnan(zmap))) == 0
        disp('Deleting files for this date and satellite because the DEM is total garbage! Moving on... ')
        delete_file = ['delete *',tifs(p).name(1:end-7),'*']; eval(delete_file);
    else
        %plot the image
        figIM = figure; set(gcf,'position',[50 650 1200 500]);
        imagesc(Y.x,Y.y,Y.z); axis xy equal; colormap gray; hold on;
        plot(S.X,S.Y,'-k','linewidth',3); hold on;
        set(gca,'xlim',[min(S.X) max(S.X)],'ylim',[min(S.Y) max(S.Y)]);
        drawnow;

        %plot the DEM
        figDEM = figure; set(gcf,'position',[50 50 1200 500]);
        imagesc(Z.x,Z.y,Z.z.raw); axis xy equal; hold on;
        colormap(gca,elev_cmap); set(gca,'clim',[nanmedian(min(zmap))-3*mad(min(zmap),1) nanmedian(min(zmap))-3*mad(min(zmap),1)+50]); cbar = colorbar;
        plot(S.X,S.Y,'-k','linewidth',3); hold on;
        set(gca,'xlim',[min(S.X) max(S.X)],'ylim',[min(S.Y) max(S.Y)]);
        drawnow;

        %check that the DEM covers a few kilometers of the melange, otherwise
        %delete the melange DEM so it isn't called again in the processing pipeline
        answer = questdlg('Does the DEM include a good cloud-free chunk of the fjord?',...
            'DEM coverage','1) Yes!','2) No!','1) Yes!');
        switch answer
            case '1) Yes!'
                elevs = Z.z.raw;

                %remind the user what to do
                disp('Iteratively generate masks to eliminate areas with obvious elevation blunders &')

                %iterative blunder removal
                q=1;
                while q
                    if q == 1
                        blunder_question = questdlg('Create blunder masks in DEM?',...
                            'Blunder ID','1) Yes!','2) No!','1) Yes!');
                        spurious_vals = zeros(size(Z.z.raw));
                    else
                        blunder_question = questdlg('Are additional masks necessary?',...
                            'Blunder ID','1) Yes!','2) No!','1) Yes!');
                    end
                    %execute blunder removal based on question response
                    switch blunder_question
                        case '1) Yes!'
                            figure(figDEM);
                            disp('Click on UL & LR corners of a box bounding the anomalous elevations to zoom in'); % Upper left, lower right.
                            [a] = ginput(2);
                            set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
                            drawnow;
                            disp('Draw mask');
                            anom_zmask = roipoly;
                            anom_zmask = double(~anom_zmask);
                            elevs = anom_zmask.*elevs; elevs(elevs==0) = NaN;
                            imagesc(Z.x,Z.y,elevs); hold on; axis xy equal;
                            colormap(gca,elev_cmap); set(gca,'clim',[nanmedian(min(zmap))-3*mad(min(zmap),1) nanmedian(min(zmap))-3*mad(min(zmap),1)+50]); cbar = colorbar;
                            plot(S.X,S.Y,'-k','linewidth',3); hold on;
                            drawnow;
                            blunders = ~anom_zmask;
                            spurious_vals = spurious_vals + blunders;
                            q = q+1;
                            set(gca,'xlim',[min(S.X) max(S.X)],'ylim',[min(S.Y) max(S.Y)]);
                            drawnow;
                        case '2) No!'
                            spurious_vals(spurious_vals>0) = 1;
                            clear blunders;
                            break
                    end
                    clear blunder_question;
                end
                Z.blunder_mask = single(spurious_vals);
                zmap(Z.blunder_mask == 1) = NaN;

                %check for even more subtle blunders
                q=1;
                while q
                    figure(figDEM); colormap(gca,elev_cmap); set(gca,'clim',[nanmedian(min(zmap))-3*mad(min(zmap),1) nanmedian(min(zmap))-3*mad(min(zmap),1)+10]); colorbar; drawnow;
                    if q == 1
                        mask_question = questdlg('Zoomed in on elevations. Any more wonky spots in DEM?',...
                            'Finer DEM Masking','1) Yes!','2) No!','1) Yes!');
                        spurious_vals = zeros(size(Z.z.raw));
                    else
                        mask_question = questdlg('Are additional masks necessary?',...
                            'Finer DEM Masking','1) Yes!','2) No!','1) Yes!');
                    end
                    %execute ocean masking based on question response
                    switch mask_question
                        case '1) Yes!'
                            figure(figDEM);
                            disp('Click on UL & LR corners of a box bounding the open ocean to zoom in');
                            [a] = ginput(2);
                            set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
                            drawnow;
                            disp('Draw mask');
                            anom_zmask = roipoly;
                            anom_zmask = double(~anom_zmask);
                            elevs = anom_zmask.*elevs; elevs(elevs==0) = NaN;
                            imagesc(Z.x,Z.y,elevs); hold on; axis xy equal;
                            colormap(gca,elev_cmap); set(gca,'clim',[nanmedian(min(zmap))-3*mad(min(zmap),1) nanmedian(min(zmap))-3*mad(min(zmap),1)+10]); cbar = colorbar;
                            plot(S.X,S.Y,'-k','linewidth',3); hold on;
                            drawnow;
                            blunders = ~anom_zmask;
                            spurious_vals = spurious_vals + blunders;
                            q = q+1;
                            set(gca,'xlim',[min(S.X) max(S.X)],'ylim',[min(S.Y) max(S.Y)]);
                            drawnow;
                        case '2) No!'
                            spurious_vals(spurious_vals>0) = 1;
                            clear blunders;
                            break
                    end
                    clear mask_question;
                end
                Z.blunder_mask = Z.blunder_mask + single(spurious_vals);
                Z.blunder_mask(Z.blunder_mask>1) = 1;

                %             %readjust elevation color scaling
                %             imagesc(double(Z.x),double(Z.y),double(Z.z.raw));
                %             colormap(gca,elev_cmap); set(gca,'clim',[nanmedian(min(Z.z.raw))-4*mad(min(Z.z.raw),1) nanmedian(min(Z.z.raw))+50]); cbar = colorbar;
                %
                %             %trace the terminus, making sure to intersect the edges of the melange mask
                %             disp('trace a line at the bottom of the terminus cliff, placing one point on each end outside of the melange outline');
                %             disp('... if the terminus isn''t visible, draw a straight line across the inland limit of observations');
                %             [term_x,term_y,~] = improfile;
                %             [Z.term.x,Z.term.y] = poly2cw(term_x,term_y);

                %save the DEM
                Z.x = single(Z.x); Z.y = single(Z.y); Z.z.raw = single(Z.z.raw);
                clear elevs;
                save(DEM_name,'Z','-v7.3'); %save as a matfile
                disp('DEM resaved as a matfile. Moving on...');
            case '2) No!'
                disp('Deleting files for this date and satellite! Moving on... ')
                delete_file = ['delete *',tifs(p).name(1:end-7),'*']; eval(delete_file);
        end
        clear answer Z Y zmap;
        close(figDEM); close(figIM); drawnow;
    end
end
disp('Done with DEM filtering & matfile creation');

%% spit out DEM dates
%find files
cd(root_dir);
mats = dir([glacier_abbrev,'*fjord-DEM.mat']); %melange DEMs provided by PGC

%pull dates
for p = 1:length(mats); DEMmat_dates(p,:) = mats(p).name(5:12);end
disp(DEMmat_dates);

%% reproject good files

%find the good DEMs that were converted to matfiles
cd(root_dir);
mats = dir('*_fjord-DEM.mat');

%pull dates
for p = 1:length(mats); DEMmat_dates(p,:) = mats(p).name(length(glacier_abbrev)+2:(length(glacier_abbrev)+2)+7);end

%calculate geoid heights & reproject
disp('Converting to orthometric heights... takes a LONG time');
for p = 1:length(mats)
    %load the DEMOCt
    disp(DEMmat_dates(p,:));
    DEM_name = [glacier_abbrev,'-',DEMmat_dates(p,:),'_fjord-DEM.mat'];
    load(DEM_name);
    if ~isfield(Z.z,'ortho')
            disp(['reformating DEM mat #',num2str(p),' of ',num2str(length(mats))]);

        %calculate geoid heights to convert from ellipsoidal heights to orthometric elevations
        disp('calculating the geoid heights for the DEM...')
        warning off;
        geoid_spacing = 100;
        [ZXgrid,ZYgrid] = meshgrid(Z.x,Z.y);
        for j = 1:geoid_spacing:size(ZXgrid,1)
            for i = 1:geoid_spacing:size(ZXgrid,2)
                [Zlon(ceil(j/geoid_spacing),ceil(i/geoid_spacing)),Zlat(ceil(j/geoid_spacing),ceil(i/geoid_spacing))] = ps2wgs(ZXgrid(j,i),ZYgrid(j,i));
                G(ceil(j/geoid_spacing),ceil(i/geoid_spacing)) = geoidheight(Zlat(ceil(j/geoid_spacing),ceil(i/geoid_spacing)),Zlon(ceil(j/geoid_spacing),ceil(i/geoid_spacing)));
            end
        end
        Z.z.geoid = single(interp2(ZXgrid(1:geoid_spacing:size(ZXgrid,1),1:geoid_spacing:size(ZXgrid,2)),ZYgrid(1:geoid_spacing:size(ZXgrid,1),1:geoid_spacing:size(ZXgrid,2)),G,ZXgrid,ZYgrid,'linear'));
        disp('converting from ellipsoidal heights to orthometric elevations...');
        Z.z.ortho = Z.z.raw - Z.z.geoid; Z.z = rmfield(Z.z,'raw');
        clear Zlat Zlon G;

        %save the file if it overlaps the melange at the terminus
        save(DEM_name,'Z','-v7.3'); %save as a matfile
        clear Z*grid Z;
        close all; drawnow;
    else
        disp('Already reformatted to orthometric heights')
    end
end