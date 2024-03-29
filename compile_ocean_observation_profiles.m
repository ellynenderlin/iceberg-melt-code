%compile ocean temperature and salinity profiles
clearvars; close all;
addpath('/users/ellynenderlin/Research/miscellaneous/general-code','/users/ellynenderlin/Research/miscellaneous/general-code/cmocean');

%% add the reference map
%load Antarctic image to plot as background for a map
Antarctic_map = 'LIMA';
cd(['/Users/ellynenderlin/Research/miscellaneous/',Antarctic_map]);
if strcmp('LIMA',Antarctic_map)
    [A,S] = readgeoraster('00000-20080319-092059124.tif');
    IM.x = S.XWorldLimits(1):S.SampleSpacingInWorldX:S.XWorldLimits(2);
    IM.y = S.YWorldLimits(2):-S.SampleSpacingInWorldY:S.YWorldLimits(1);
    IM.z=single(A); IM.z = IM.z./255;
else %assume REMA
    IM.x =S.XWorldLimits(1)+0.5*S.CellExtentInWorldX:S.CellExtentInWorldX:S.XWorldLimits(2)-0.5*S.CellExtentInWorldX;
    IM.y = S.YWorldLimits(2)-0.5*S.CellExtentInWorldY:-S.CellExtentInWorldY:S.YWorldLimits(1)+0.5*S.CellExtentInWorldY;
    IM.z=single(A);
end
clear A S;


%% compile CTD & APB data (only run once to create Antarctic-ocean-data.mat)

%navigate to the CTD and APB data
CTD_path = '/users/ellynenderlin/Research/NSF_Antarctic-Icebergs/CTD_Antarctica/';
cd(CTD_path);

%locate the data files
CTD_files = dir('*CTD.mat'); APB_files = dir('*APB.mat'); ARGO_files = dir('*ARGO.mat');

%load the CTD data
overviewfig = figure; set(overviewfig,'position',[50 50 1200 1200]);
cmap = colormap(jet(length(CTD_files)+length(APB_files)+length(ARGO_files)));
sub1 = subplot(1,2,1); sub2 = subplot(1,2,2);
set(sub2,'position',[0.7 0.11 0.25 0.85]); set(sub1,'position',[0.11 0.11 0.55 0.85]); 
subplot(sub1); imagesc(IM.x,IM.y,IM.z); colormap gray; hold on; axis xy equal
colormap(gray(10001));
for i = 1:length(CTD_files)
    load_file = ['load ',CTD_files(i).name]; eval(load_file); %info_sub describes variables
    
    %extract region name
    refs = strfind(CTD_files(i).name,'_');
    if length(refs) > 1
        region_name(i) = {CTD_files(i).name(refs(1)+1:refs(2)-1)};
    else
        region_name(i) = {CTD_files(i).name(1:refs(1)-1)};
    end
    CTD(i).name = string(region_name(i));
    
    if ~strcmp(CTD_files(i).name,'WesternPeninsula_CTD.mat')
        %add positional data to structure
        CTD(i).lat = latitude; %latitude (degrees N)
        CTD(i).lon = longitude; %longitude (degrees E)
        [CTD(i).x,CTD(i).y] = wgs2ps(CTD(i).lon,CTD(i).lat,'StandardParallel',-71,'StandardMeridian',0);
        CTD(i).depth = depth; %depth (m)
        
        %add dates to structure
        CTD(i).time = date_datenum; %matlab dates
        [YYYY,MM,DD,~,~,~] = datevec(CTD(i).time); decidate = NaN(length(YYYY),1); 
        %convert years to strings
        years = num2str(YYYY);
        %convert months to strings
        for j = 1:length(MM); if MM(j) < 10; mos(j,:)= ['0',num2str(MM(j))]; else mos(j,:) = num2str(MM(j)); end; end
        %convert days to strings
        for j = 1:length(DD); if DD(j) < 10; days(j,:)= ['0',num2str(DD(j))]; else days(j,:) = num2str(DD(j)); end; end
        datestrings = [years,mos,days];
        %convert to decimal date
        for j = 1:length(years); decidate(j) = convert_to_decimaldate(datestrings(j,:)); end
        CTD(i).date = decidate; %decimal date
        clear YYYY MM DD decidate datestrings years mos days;
        
        %add data to structure
        CTD(i).P = []; %pressure (dbar)
        CTD(i).T = temperature; %in situ or conservative temp? (C)
        CTD(i).sal = salinity; %practical or absolute salinity?
    else %W Peninsula is a combination of different datasets (including CTD and APB)
        %add positional data to structure
        CTD(i).lat = lat_sub; %latitude (degrees N)
        CTD(i).lon = lon_sub+360; %longitude (degrees E)
        [CTD(i).x,CTD(i).y] = wgs2ps(CTD(i).lon,CTD(i).lat,'StandardParallel',-71,'StandardMeridian',0);
        CTD(i).depth = repmat(z_sub,1,length(time_sub)); %height (m)
        
        %add dates to structure
        CTD(i).time = time_sub; %matlab dates
        [YYYY,MM,DD,~,~,~] = datevec(CTD(i).time); decidate = NaN(length(YYYY),1); 
        %convert years to strings
        years = num2str(YYYY');
        %convert months to strings
        for j = 1:length(MM); if MM(j) < 10; mos(j,:)= ['0',num2str(MM(j))]; else mos(j,:) = num2str(MM(j)); end; end
        %convert days to strings
        for j = 1:length(DD); if DD(j) < 10; days(j,:)= ['0',num2str(DD(j))]; else days(j,:) = num2str(DD(j)); end; end
        datestrings = [years,mos,days];
        %convert to decimal date
        for j = 1:length(years); decidate(j) = convert_to_decimaldate(datestrings(j,:)); end
        CTD(i).date = decidate; %decimal date
        clear YYYY MM DD decidate datestrings years mos days;
        
        %add data to structure
        CTD(i).P = p_sub; %pressure (dbar)
        CTD(i).T = temp_sub; %in situ or conservative temp? (C)
        CTD(i).sal = salt_sub; %practical or absolute salinity?
    end
    subplot(sub1);
    plot(CTD(i).x,CTD(i).y,'x','color',cmap(i,:)); hold on; drawnow;
    subplot(sub2);
    plot(CTD(i).T,CTD(i).depth,'x','color',cmap(i,:)); hold on; drawnow;
    pl(i) = plot(CTD(i).T(1),CTD(i).depth(1),'x','color',cmap(i,:)); hold on; drawnow;
    
    clearvars -except CTD* APB* ARGO* cmap sub* pl region_name *days*;
end

%load the APB data
for i = 1:length(APB_files)
    load_file = ['load ',APB_files(i).name]; eval(load_file); %info_sub describes variables
    
    %extract region name
    refs = strfind(APB_files(i).name,'_');
    if length(refs) > 1
        region_name(i+length(CTD_files)) = {APB_files(i).name(refs(1)+1:refs(2)-1)};
    else
        region_name(i+length(CTD_files)) = {APB_files(i).name(1:refs(1)-1)};
    end
    APB(i).name = string(region_name(i+length(CTD_files)));
    
    %add positional data to structure
    APB(i).lat = latitude; %latitude (degrees N)
    APB(i).lon = longitude; %longitude (degrees E)
    [APB(i).x,APB(i).y] = wgs2ps(APB(i).lon,APB(i).lat,'StandardParallel',-71,'StandardMeridian',0);
    APB(i).depth = depth; %depth (m)
    
    %add dates to structure
    APB(i).time = date_datenum; %matlab dates
    [YYYY,MM,DD,~,~,~] = datevec(APB(i).time); decidate = NaN(length(YYYY),1);
    %convert years to strings
    years = num2str(YYYY);
    %convert months to strings
    for j = 1:length(MM); if MM(j) < 10; mos(j,:)= ['0',num2str(MM(j))]; else mos(j,:) = num2str(MM(j)); end; end
    %convert days to strings
    for j = 1:length(DD); if DD(j) < 10; days(j,:)= ['0',num2str(DD(j))]; else days(j,:) = num2str(DD(j)); end; end
    datestrings = [years,mos,days];
    %convert to decimal date
    for j = 1:length(years); decidate(j) = convert_to_decimaldate(datestrings(j,:)); end
    APB(i).date = decidate; %decimal date
    clear YYYY MM DD decidate datestrings years mos days;
    
    %add data to structure
    APB(i).P = []; %pressure (dbar)
    APB(i).T = temperature; %in situ or conservative temp? (C)
    for j = 1:size(APB(i).T,2)
       if max(APB(i).T(:,j)) > 5;  APB(i).T(:,j) = NaN; end %remove strangely high temps when I assume the seals are beached
    end
    APB(i).sal = salinity; %practical or absolute salinity?
    subplot(sub1);
    plot(APB(i).x,APB(i).y,'x','color',cmap(i+length(CTD_files),:)); hold on; drawnow;
    subplot(sub2);
    plot(APB(i).T,APB(i).depth,'x','color',cmap(i+length(CTD_files),:)); hold on; drawnow;
    pl(i+length(CTD_files)) = plot(APB(i).T(1),APB(i).depth(1),'x','color',cmap(i+length(CTD_files),:)); hold on; drawnow;
    
    clearvars -except CTD* APB* ARGO* cmap sub* pl region_name *days*;
end

%load the ARGO data
for i = 1:length(ARGO_files)
    load_file = ['load ',ARGO_files(i).name]; eval(load_file); %info_sub describes variables
    
    %extract region name
    refs = strfind(ARGO_files(i).name,'_');
    if length(refs) > 1
        region_name(i+length(CTD_files)+length(APB_files)) = {ARGO_files(i).name(refs(1)+1:refs(2)-1)};
    else
        region_name(i+length(CTD_files)+length(APB_files)) = {ARGO_files(i).name(1:refs(1)-1)};
    end
    ARGO(i).name = string(region_name(i+length(CTD_files)+length(APB_files)));
    
    %add positional data to structure
    ARGO(i).lat = latitude; %latitude (degrees N)
    ARGO(i).lon = longitude; %longitude (degrees E)
    [ARGO(i).x,ARGO(i).y] = wgs2ps(ARGO(i).lon,ARGO(i).lat,'StandardParallel',-71,'StandardMeridian',0);
    ARGO(i).depth = pressure; %depth (m)... inferred directly from ARGO depth (dbar)
    
    %add dates to structure
    ARGO(i).time = date_datenum; %matlab dates
    [YYYY,MM,DD,~,~,~] = datevec(ARGO(i).time); decidate = NaN(length(YYYY),1);
    %convert years to strings
    years = num2str(YYYY);
    %convert months to strings
    for j = 1:length(MM); if MM(j) < 10; mos(j,:)= ['0',num2str(MM(j))]; else mos(j,:) = num2str(MM(j)); end; end
    %convert days to strings
    for j = 1:length(DD); if DD(j) < 10; days(j,:)= ['0',num2str(DD(j))]; else days(j,:) = num2str(DD(j)); end; end
    datestrings = [years,mos,days];
    %convert to decimal date
    for j = 1:length(years); decidate(j) = convert_to_decimaldate(datestrings(j,:)); end
    ARGO(i).date = decidate; %decimal date
    clear YYYY MM DD decidate datestrings years mos days;
    
    %add data to structure
    ARGO(i).P = []; %pressure (dbar)
    ARGO(i).T = temperature; %in situ or conservative temp? (C)
    ARGO(i).sal = salinity; %practical or absolute salinity?
    subplot(sub1);
    plot(ARGO(i).x,ARGO(i).y,'x','color',cmap(i+length(CTD_files)+length(APB_files),:)); hold on; drawnow;
    subplot(sub2);
    plot(ARGO(i).T,ARGO(i).depth,'x','color',cmap(i+length(CTD_files)+length(APB_files),:)); hold on; drawnow;
    pl(i+length(CTD_files)+length(APB_files)) = plot(ARGO(i).T(1),ARGO(i).depth(1),'x','color',cmap(i+length(CTD_files)+length(APB_files),:)); hold on; drawnow;
    
    clearvars -except CTD* APB* ARGO* cmap sub* pl region_name *days*;
end
%plot
subplot(sub2);
leg = legend(pl,region_name);
subplot(sub1);
set(gca,'xlim',[-28e5 28e5],'xtick',[-32e5:8e5:32e5],'xticklabel',[-3200:800:3200],...
    'ylim',[-24e5 24e5],'ytick',[-24e5:8e5:24e5],'yticklabel',[-2400:800:2400]); grid off;


%save the data
save([CTD_path,'Antarctic-ocean-data.mat'],'CTD','APB','ARGO','-v7.3');
disp('All available ocean temperature data saved to one Matlab file');
