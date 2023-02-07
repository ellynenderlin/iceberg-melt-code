function [SMB,temp,smb] = read_RACMO_Greenland_netcdfs(DEM1_date,DEM2_date,dir_RACMO)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


%convert input DEM dates from strings to decimal dates
startyr = convert_to_decimaldate(DEM1_date);
endyr = convert_to_decimaldate(DEM2_date);

%find all the RACMO files
RACMOs = dir([dir_RACMO,'smb*.nc']);
for j = 1:length(RACMOs)
    %grab year from filename and convert to number
    RACMOyr(j,1) = str2num(RACMOs(j).name(9:12)); %file name indexing for date depends on the variable!
    
    %use month indexing in filename to augment the year so that it is
    %decimal year based on the first date of the file
    if contains(RACMOs(j).name,'AMJ')
        RACMOyr(j,1) = RACMOyr(j,1) + 0.25;
    elseif contains(RACMOs(j).name,'JAS')
        RACMOyr(j,1) = RACMOyr(j,1) + 0.5;
    elseif contains(RACMOs(j).name,'OND')
        RACMOyr(j,1) = RACMOyr(j,1) + 0.75;
    end
end

%sort RACMO data by date
[sortedyr,sortinds] = sort(RACMOyr);

%identify the files within the date range
startind = find(sortedyr <= startyr,1,'last');
endind = find(sortedyr <= endyr,1,'last');
inds = startind:1:endind;

%load files in chronological order & compile data
for j = 1:length(inds)
    %identify appropriate RACMO grid cell using data from the first file
    if j == 1
%         %load lat,lon
%         smb_lat = ncread([dir_RACMO,RACMOs(sortinds(inds(j))).name],'LAT');
%         smb_lon = ncread([dir_RACMO,RACMOs(sortinds(inds(j))).name],'LON');

        %IMPORTANT NOTE: RACMO2.3_Antarctica reference grid cell
        %identification relied on the reprojection of lat,lon to Polar
        %Stereo coordinates. As a result, when pulling the data from RACMO,
        %the y,x references were used rather than the x,y like below!
        
        %Polar Stereo X,Y
        smb_x = ncread([dir_RACMO,RACMOs(sortinds(inds(j))).name],'x');
        smb_y = ncread([dir_RACMO,RACMOs(sortinds(inds(j))).name],'y');
        [xgrid,ygrid] = meshgrid(smb_x,smb_y); clear smb_x smb_y;
        smb_x = xgrid; smb_y = ygrid; clear xgrid ygrid;
        
        %grab the variables
        smb = ncread([dir_RACMO,RACMOs(sortinds(inds(j))).name],'smb_rec'); smb(smb==-9999)=NaN;
        %         runoff = ncread([dir_RACMO,'SMB',RACMOs(sortinds(inds(j))).name(7:end)],'runoffcorr'); smb(smb==-9999)=NaN;
        %         airtemp = ncread([dir_RACMO,'T2m',RACMOs(sortinds(inds(j))).name(7:end)],'t2m'); airtemp(airtemp==-9999)=NaN;
        %         icetemp = ncread([dir_RACMO,'T2m',RACMOs(sortinds(inds(j))).name(7:end)],'t2m'); icetemp(icetemp==-9999)=NaN; %assume the ice temp matches the avg long-term temp (for creep estimates)
        smb_time = ncread([dir_RACMO,RACMOs(sortinds(inds(j))).name],'time');
        
        %calculate the x&y distance between the target of interest (berg lon, berg lat) and each RACMO grid cell
        x_diff = abs(berg_x - smb_x);
        y_diff = abs(berg_y - smb_y);
        diff_map = sqrt(y_diff.^2+x_diff.^2); diff_map(isnan(squeeze(nanmean(smb(:,:,:),3)))) = NaN;%solve for the distance vector using the x&y distances
        RACMO_ref = find(diff_map==min(min(diff_map))); %find the minimum distance (reference for your grid cell is output)
        [RACMOy RACMOx] = ind2sub(size(diff_map),RACMO_ref); %convert cell reference to an x- and y-cell index
        disp(['RACMO x-reference = ',num2str(RACMOx),' & y-reference = ',num2str(RACMOy)]);
        
        %adjust RACMO reference grid cell if necessary
        disp('Adjust coordinates (if necessary) to extract surface melt estimates');
        figure; set(gcf,'position',[100 100 700 700]);
        smb_cmap = colormap(jet(10001)); smb_cmap(1,:) = [1 1 1];
        smb_map = max(smb,[],3).*86500./1000; smb_map(isnan(smb_map)) = min(smb_map(~isnan(smb_map)))-1;
        
        %plot map with actual coordinates to confirm proper referencing
%         imagesc(smb_x(1,:),smb_y(:,1),smb_map'); colormap(gca,smb_cmap); hold on; axis xy equal;
%         plot(berg_x,berg_y,'xc','linewidth',2,'markersize',20);
%         plot(smb_x(1,RACMOx),smb_y(RACMOy,1),'+m','linewidth',2,'markerfacecolor','none','markersize',20);
%         set(gca,'xlim',[smb_x(1,RACMOx)-50000 smb_x(1,RACMOx)+50000],...
%             'ylim',[smb_y(RACMOy,1)-50000 smb_y(RACMOy,1)+50000]);

        %plot map with grid indices to identify nearest neighboring
        %reference grid cell with data
        imagesc(smb_map'); colormap(gca,smb_cmap); hold on; axis xy equal;
        plot(RACMOx,RACMOy,'+m','linewidth',2,'markerfacecolor','none','markersize',20);
        disp(['smb at RACMO pixel = ',num2str(max(smb(RACMOx,RACMOy,:)).*86500./1000),' mm w.e.']);
        set(gca,'clim',[0 max(smb(RACMOx,RACMOy,:)).*86500./1000+2*std(smb((smb~=0))).*86500./1000]); cbar = colorbar; 
        set(get(cbar,'ylabel'),'string','SMB (m w.e. per day');
        set(gca,'xlim',[RACMOx-20 RACMOx+20],...
            'ylim',[RACMOy-20 RACMOy+20]);
        
        %adjust coordinates as necessary
        disp('Modify coordinates if the marker is in a region with no data');
        adjustment = questdlg('Do the coordinates need to be modified?',...
            'RACMO coordinate adjustment','yes','no','yes');
        switch adjustment
            case 'yes'
                disp('To change coordinates, identify new coordinates using tick marks and type RACMOx=XX; RACMOy=YY; dbcont ');
                keyboard
            case 'no'
                disp('Using default coordinates');
        end
        clear adjustment;
        close all; drawnow;
    else
        smbsub = ncread([dir_RACMO,RACMOs(sortinds(inds(j))).name],'smb_rec'); smbsub(smbsub==-9999)=NaN;
        %         runoffsub = ncread([dir_RACMO,'SMB',RACMOs(sortinds(inds(j))).name(7:end)],'runoffcorr'); runoffsub(runoffsub==-9999)=NaN;
        %         airtempsub = ncread([dir_RACMO,'T2m',RACMOs(sortinds(inds(j))).name(7:end)],'t2m'); airtempsub(airtempsub==-9999)=NaN;
        %         icetempsub = ncread([dir_RACMO,'T2m',RACMOs(sortinds(inds(j))).name(7:end)],'t2m'); icetempsub(icetempsub==-9999)=NaN; %assume the ice temp matches the avg long-term temp (for creep estimates)
        smbsub_time = ncread([dir_RACMO,RACMOs(sortinds(inds(j))).name],'time');
        
        %concatenate the data to create continuous timeseries
        smb = cat(3,smb,smbsub);
%         runoff = cat(3,runoff,runoffsub);
%         airtemp = cat(3,airtemp,airtempsub);
%         icetempsub = cat(3,icetempsub,icetempsub);
        
        %PICK UP HERE... NEED THESE TO BE DATES SINCE SOME STANDARD TIME,
        %NOT 0-90 WITH RESTARTS FOR EACH FILE
        smb_time = cat(3,smb_time,smbsub_time);
    
    end
    
    
    
    
    
    
end




end