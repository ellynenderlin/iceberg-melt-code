function [dt,berg_T,berg_runoff] = extract_MAR_params(dir_SMB,geography,berg_x,berg_y,berg_dates)
%%%extract SMB parameters from Greenland MAR
if geography ~=0 
    error('Use RACMO data for Antarctica!');
end

%convert input DEM dates from strings to decimal dates
startyr = convert_to_decimaldate(berg_dates(1,:));
endyr = convert_to_decimaldate(berg_dates(2,:));

%find all the annual MAR files
MARs = dir([dir_SMB,'MAR*.nc']);
for j = 1:length(MARs)
    MARyrs(j) = str2num(MARs(j).name(end-6:end-3));
end

%identify the files within the date range
startind = find(MARyrs <= startyr,1,'last');
endind = find(MARyrs <= endyr,1,'last');
inds = startind:1:endind; clear startind endind;

%loop through annual MAR files within the date range & concatenate
smb_stack = []; runoff_stack = []; airtemp_stack = []; icetemp_stack = []; date_stack = [];
for j = inds
    
    %create time stamp vector
    smb_yr = ncread([dir_SMB,MARs(j).name],'YYYY');
    smb_mo = ncread([dir_SMB,MARs(j).name],'MM');
    smb_da = ncread([dir_SMB,MARs(j).name],'DD');
    for k = 1:length(smb_yr)
        if smb_mo(k)<10
            if smb_da(k) < 10
                [smb_time(k)] = convert_to_decimaldate([num2str(smb_yr(k,:)),'0',num2str(smb_mo(k,:)),'0',num2str(smb_da(k,:))]);
            else
                [smb_time(k)] = convert_to_decimaldate([num2str(smb_yr(k,:)),'0',num2str(smb_mo(k,:)),num2str(smb_da(k,:))]);
            end
        else
            if smb_da(k) < 10
                [smb_time(k)] = convert_to_decimaldate([num2str(smb_yr(k,:)),num2str(smb_mo(k,:)),'0',num2str(smb_da(k,:))]);
            else
                [smb_time(k)] = convert_to_decimaldate([num2str(smb_yr(k,:)),num2str(smb_mo(k,:)),num2str(smb_da(k,:))]);
            end
        end
    end
    
    %grab the variables
    GRISmask = squeeze(ncread([dir_SMB,MARs(j).name],'MSK')); 
    smb = squeeze(ncread([dir_SMB,MARs(j).name],'SMB')); smb(smb==-1.0000e+19) = NaN; % surface mass balance (mmWE/day)
    runoff = squeeze(ncread([dir_SMB,MARs(j).name],'RU')); runoff(runoff==-1.0000e+19) = NaN; %meltwater+rain runoff (mmWE/day)
    airtemp = squeeze(ncread([dir_SMB,MARs(j).name],'TT')); airtemp(airtemp==-1.0000e+19) = NaN; %surface temp (degrees C)
    %mask-out water
    GRISmask(GRISmask==0) = NaN; GRISmask(~isnan(GRISmask)) = 1;
    for k = 1:size(smb,3)
        smb(:,:,k) = smb(:,:,k).*GRISmask;
        runoff(:,:,k) = runoff(:,:,k).*GRISmask;
        airtemp(:,:,k) = airtemp(:,:,k).*GRISmask;
    end
    icetemp = airtemp;
    
    %extract coordinates & identify nearest MAR grid cell during first loop
    if j == inds(1)
        %Polar Stereo X,Y
        smb_x = ncread([dir_SMB,MARs(j).name],'x'); smb_x = 1000*smb_x; %convert from km to m
        smb_y = ncread([dir_SMB,MARs(j).name],'y'); smb_y = 1000*smb_y; %convert from km to m
        [xgrid,ygrid] = meshgrid(smb_x,smb_y); %clear smb_x smb_y;
        %     smb_x = xgrid; smb_y = ygrid; clear xgrid ygrid;
        
        %calculate the x&y distance between the target of interest (berg lon, berg lat) and each MAR grid cell
        x_diff = abs(repmat(berg_x,size(xgrid)) - xgrid);
        y_diff = abs(repmat(berg_y,size(ygrid)) - ygrid);
        diff_map = sqrt(y_diff.^2+x_diff.^2); diff_map(isnan(squeeze(nanmean(smb(:,:,:),3)))) = NaN;%solve for the distance vector using the x&y distances
        MAR_ref = find(diff_map==min(min(diff_map))); %find the minimum distance (reference for your grid cell is output)
        [MARy MARx] = ind2sub(size(diff_map),MAR_ref); %convert cell reference to an x- and y-cell index
        disp(['MAR x-reference = ',num2str(MARx),' & y-reference = ',num2str(MARy)]);
        
        %adjust MAR reference grid cell if necessary
        disp('Adjust coordinates (if necessary) to extract surface melt estimates');
        figure; set(gcf,'position',[100 100 700 700]);
        smb_cmap = cmocean('thermal',10001); smb_cmap(1,:) = [1 1 1];
        smb_map = sum(smb,3); 
        disp(['annual smb at MAR pixel = ',num2str(smb_map(MARx,MARy)./1000),' m w.e.']);
        
        %plot map with grid indices to identify nearest neighboring
        %reference grid cell with data
%         smb_map(isnan(smb_map)) = min(smb_map(~isnan(smb_map)))-1; %replace NaN with smallest observed value
        imagesc(smb_map'); colormap(gca,smb_cmap); hold on; axis xy equal;
        plot(MARx,MARy,'+k','linewidth',2,'markerfacecolor','none','markersize',20);
%         set(gca,'clim',[smb_map(MARx,MARy)-1*nanstd(smb_map(~isnan(smb_map))) smb_map(MARx,MARy)+1*nanstd(smb_map(~isnan(smb_map)))]);
        cbar = colorbar; set(get(cbar,'ylabel'),'string','SMB (mm w.e. per day');
        set(gca,'xlim',[MARx-20 MARx+20],...
            'ylim',[MARy-20 MARy+20]);

        %adjust coordinates as necessary
        disp('Modify coordinates if the marker is in a region with no data');
        adjustment = questdlg('Do the coordinates need to be modified?',...
            'MAR coordinate adjustment','yes','no','yes');
        switch adjustment
            case 'yes'
                disp('To change coordinates, identify new coordinates using tick marks and type MARx=XX; MARy=YY; dbcont ');
                keyboard
            case 'no'
                disp('Using default coordinates');
        end
        clear adjustment;
        close all; drawnow;
    end
    
    %identify the date range to use for each annual file
    startind = find(smb_time >= startyr,1,'first');
    endind = find(smb_time <= endyr,1,'last');
    
    %concatenate data from the different annual files
    smb_stack = cat(3,smb_stack,smb(MARx,MARy,startind:endind)); 
    runoff_stack = cat(3,runoff_stack,runoff(MARx,MARy,startind:endind)); 
    airtemp_stack = cat(3,airtemp_stack,airtemp(MARx,MARy,startind:endind)); 
    icetemp_stack = cat(3,icetemp_stack,icetemp(MARx,MARy,:));
    date_stack = cat(1,date_stack,smb_time(startind:endind)');
    
    %clear annual variables
    clear smb_yr smb_mo smb_da smb_time startind endind GRISmask smb runoff airtemp icetemp;
end
smb = squeeze(smb_stack); runoff = squeeze(runoff_stack); airtemp = squeeze(airtemp_stack); icetemp = squeeze(icetemp_stack); MARdates = date_stack;
clear *_stack;


%calculate the time separation between DEMs in terms of
%decimal years (ddays) & decimal days (days)
to = berg_dates(1,:); tf = berg_dates(2,:);
dt = datenum(tf(1:12),'yyyymmddHHMM') - datenum(to(1:12),'yyyymmddHHMM');
days = ones(1,(datenum(tf(1:8),'yyyymmdd') - datenum(to(1:8),'yyyymmdd'))+1); days(2:end-1) = 1; 
days(1) = ceil(datenum(to,'yyyymmddHHMMSS'))-datenum(to,'yyyymmddHHMMSS'); 
days(end) = datenum(tf,'yyyymmddHHMMSS')-floor(datenum(tf,'yyyymmddHHMMSS'));

%estimate surface melting using MAR runoff (mm w.e.)
berg_runoff = nansum(days'.*runoff)/1000; %surface meltwater that runs off (convert from mm w.e. per day to total m w.e.)
%estimate the ice temperature as the average annual air temperature (doesn't account for advection of colder ice and melt/refreezing at the surface and/or submarine interface)
berg_T = 273+nanmean(icetemp); % air temp (Kelvin)



end