%compile ocean temperature and salinity profiles
clearvars; close all;
addpath('/users/ellynenderlin/miscellaneous/general-code','/users/ellynenderlin/miscellaneous/general-code/cmocean');


%% compile CTD & APB data (only run once to create Antarctic-ocean-data.mat)

%navigate to the CTD and APB data
CTD_path = '/users/ellynenderlin/Research/NSF_Antarctic-Icebergs/CTD_Antarctica/';
cd(CTD_path);

%locate the data files
CTD_files = dir('*CTD.mat'); APB_files = dir('*APB.mat'); 

%set up days of year for leap and non-leap years for date conversions
modays = [31 28 31 30 31 30 31 31 30 31 30 31]; cumdays = [0 cumsum(modays(1:end-1))];
leap_modays = [31 29 31 30 31 30 31 31 30 31 30 31]; leap_cumdays = [0 cumsum(leap_modays(1:end-1))];

%load the CTD data
overviewfig = figure; set(overviewfig,'position',[50 50 1200 1200]);
cmap = colormap(jet(length(CTD_files)+length(APB_files)));
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
        CTD(i).lat = latitude; %latitude (degrees N)
        CTD(i).lon = longitude; %longitude (degrees E)
        [CTD(i).x,CTD(i).y] = wgs2ps(CTD(i).lon,CTD(i).lat,'StandardParallel',-71,'StandardMeridian',0);
        CTD(i).depth = depth; %depth (m)
        CTD(i).time = date_datenum; %matlab dates
        [YYYY,MM,DD,~,~,~] = datevec(CTD(i).time); decidate = NaN(length(YYYY),1); 
        decidate(mod(YYYY,4)~=0) = YYYY(mod(YYYY,4)~=0) + ((cumdays(MM(mod(YYYY,4)~=0))'+DD(mod(YYYY,4)~=0))./sum(modays));
        decidate(mod(YYYY,4)==0) = YYYY(mod(YYYY,4)==0) + ((leap_cumdays(MM(mod(YYYY,4)==0))'+DD(mod(YYYY,4)==0))./sum(leap_modays));
        CTD(i).date = decidate; %decimal date
        clear YYYY MM DD hh mm ss decidate;
        CTD(i).P = []; %pressure (dbar)
        CTD(i).T = temperature; %in situ or conservative temp? (C)
        CTD(i).sal = salinity; %practical or absolute salinity?
    else %W Peninsula is a combination of different datasets (including CTD and APB)
        CTD(i).lat = lat_sub; %latitude (degrees N)
        CTD(i).lon = lon_sub+360; %longitude (degrees E)
        [CTD(i).x,CTD(i).y] = wgs2ps(CTD(i).lon,CTD(i).lat,'StandardParallel',-71,'StandardMeridian',0);
        CTD(i).depth = repmat(z_sub,1,length(time_sub)); %height (m)
        CTD(i).time = time_sub; %matlab dates
        [YYYY,MM,DD,~,~,~] = datevec(CTD(i).time); decidate = NaN(length(YYYY),1); 
        decidate(mod(YYYY,4)~=0) = YYYY(mod(YYYY,4)~=0) + ((cumdays(MM(mod(YYYY,4)~=0))+DD(mod(YYYY,4)~=0))./sum(modays));
        decidate(mod(YYYY,4)==0) = YYYY(mod(YYYY,4)==0) + ((leap_cumdays(MM(mod(YYYY,4)==0))+DD(mod(YYYY,4)==0))./sum(leap_modays));
        CTD(i).date = decidate; %decimal date
        clear YYYY MM DD hh mm ss decidate;
        CTD(i).P = p_sub; %pressure (dbar)
        CTD(i).T = temp_sub; %in situ or conservative temp? (C)
        CTD(i).sal = salt_sub; %practical or absolute salinity?
    end
    subplot(sub1);
    plot(CTD(i).x,CTD(i).y,'x','color',cmap(i,:)); hold on; drawnow;
    subplot(sub2);
    plot(CTD(i).T,CTD(i).depth,'x','color',cmap(i,:)); hold on; drawnow;
    pl(i) = plot(CTD(i).T(1),CTD(i).depth(1),'x','color',cmap(i,:)); hold on; drawnow;
    
    clearvars -except CTD* APB* cmap sub* pl region_name *days*;
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
    
    APB(i).lat = latitude; %latitude (degrees N)
    APB(i).lon = longitude; %longitude (degrees E)
    [APB(i).x,APB(i).y] = wgs2ps(APB(i).lon,APB(i).lat,'StandardParallel',-71,'StandardMeridian',0);
    APB(i).depth = depth; %depth (m)
    APB(i).time = date_datenum; %matlab dates
    [YYYY,MM,DD,~,~,~] = datevec(APB(i).time); decidate = NaN(length(YYYY),1);
    decidate(mod(YYYY,4)~=0) = YYYY(mod(YYYY,4)~=0) + ((cumdays(MM(mod(YYYY,4)~=0))'+DD(mod(YYYY,4)~=0))./sum(modays));
    decidate(mod(YYYY,4)==0) = YYYY(mod(YYYY,4)==0) + ((leap_cumdays(MM(mod(YYYY,4)==0))'+DD(mod(YYYY,4)==0))./sum(leap_modays));
    APB(i).date = decidate; %decimal date
    clear YYYY MM DD hh mm ss decidate;
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
    
    clearvars -except CTD* APB* cmap sub* pl region_name *days*;
end
subplot(sub2);
leg = legend(pl,region_name);
subplot(sub1);
set(gca,'xlim',[-28e5 28e5],'xtick',[-32e5:8e5:32e5],'xticklabel',[-3200:800:3200],...
    'ylim',[-24e5 24e5],'ytick',[-24e5:8e5:24e5],'yticklabel',[-2400:800:2400]); grid off;

% %plot time series of temperature profiles for each dataset
% for i = 1:length(CTD)
%     %identify unique years
%     yrs = unique(floor(CTD(i).date));
%     
%     %set up the figure
%     figure; yr_cmap = colormap(jet(length(yrs)));
%     for j = 1:length(CTD(i).date)
%        plot(CTD(i).T(:,j),CTD(i).depth(:,j),'-','color',yr_cmap(find(yrs==floor(CTD(i).date(j))),:)); hold on;
%     end
%     drawnow;
%     xlabel('Temperature (C)'); ylabel('Depth (m)'); title(region_name(i));
%     clear yrs yr_cmap;
% end
% for i = 1:length(APB)
%     %identify unique years
%     yrs = unique(floor(APB(i).date));
%     
%     %set up the figure
%     figure; yr_cmap = colormap(jet(length(yrs)));
%     for j = 1:length(APB(i).date)
%        plot(APB(i).T(:,j),APB(i).depth(:,j),'-','color',yr_cmap(find(yrs==floor(APB(i).date(j))),:)); hold on;
%     end
%     drawnow;
%     xlabel('Temperature (C)'); ylabel('Depth (m)'); title(region_name(i+length(CTD)));
%     clear yrs yr_cmap;
% end

%save the data
save([CTD_path,'Antarctic-ocean-data.mat'],'CTD','APB','-v7.3');