%PLOT_MELTRATE_FIGS_UPDATED: Generate maps and melt rate subplots for all
%Antarctic iceberg melt datasets. Includes ocean observation data in
%figures. Exports concatenated data for all sites to a table.


%% Set-up workspace
clearvars; close all; drawnow;
addpath('/users/ellynenderlin/Research/miscellaneous/general-code','/users/ellynenderlin/Research/miscellaneous/general-code/cmocean');

root_dir = '/Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt/';
plot_yrs = []; avgx = []; avgy = []; region = []; warning off;
rho_sw = 1026; %sea water density in kg m^-3

%make colormaps
yrs = [2013:1:2022]; %plot_color = colormap(parula(length(yrs)+1)); plot_color = plot_color(1:length(yrs),:);
% plot_color = [127,59,8;179,88,6;224,130,20;253,184,99;254,224,182;216,218,235;178,171,210;128,115,172;84,39,136;45,0,75]/255; %orange-purple colormap
plot_color = cmocean('ice',length(yrs)+1); plot_color = plot_color(1:end-1,:); %colormap for emphasizing different years
lowmelt_cmap = flipud(colormap(hot(40))); highmelt_cmap = flipud(colormap(hot(100)));
% plot_cmap = colormap(parula(15)); %colormap for emphasizing different locations
close all;

%if subplots for each region are organized 2 rows by 8 columns
region = [{'Edgeworth-LarsenA'},{'Crane-LarsenB'},{'Ronne'},{'Filchner'},{'Amery'},{'Totten'},{'Mertz'},...
    {'Thwaites'},{'Ferrigno-Eltanin'},{'Seller-Bugge'},{'Heim-Marguerite'},{'Widdowson-Biscoe'},{'Cadman-Palmer'},{'Blanchard-Danco'},{'Leonardo-Danco'}]; %site directory names
disp_names = [{'i) Edgeworth'},{'j) Crane'},{'k) Ronne'},{'l) Filchner'},{'m) Polar Times'},{'n) Totten'},{'o) Mertz'},...
    {'h) Thwaites'},{'g) Ferrigno'},{'f) Seller'},{'e) Heim'},{'d) Widdowson'},{'c) Cadman'},{'b) Blanchard'},{'a) Leonardo'}]; %geographically-organized figure labels
leg_names = [{'Edgeworth'},{'Crane'},{'Ronne'},{'Filchner'},{'Polar Times'},{'Totten'},{'Mertz'},...
    {'Thwaites'},{'Ferrigno'},{'Seller'},{'Heim'},{'Widdowson'},{'Cadman'},{'Blanchard'},{'Leonardo'}]; %generic figure labels
map_marker = 's'; %marker symbol for maps
plot_marker = 's'; %marker symbol for scatterplots
symbol_size = 4; %symbols size
plot_loc = [2,4,6,8,10,12,14,15,13,11,9,7,5,3,1]; %specifies the subplot locations for the study sites listed as "region", "disp_names", and "leg_names" above
region_colors = [77,172,38; 77,172,38; 184,225,134; 184,225,134; 184,225,134; 184,225,134; 184,225,134;...
    241,182,218; 241,182,218; 208,28,139; 208,28,139; 208,28,139; 208,28,139; 208,28,139; 208,28,139]./255;

%days of year (needed for datestring to decimal year conversion)
modays_norm = [31 28 31 30 31 30 31 31 30 31 30 31];
cumdays_norm = cumsum(modays_norm); cumdays_norm = [0 cumdays_norm(1:11)];
modays_leap = [31 29 31 30 31 30 31 31 30 31 30 31];
cumdays_leap = cumsum(modays_leap); cumdays_leap = [0 cumdays_leap(1:11)];


%% create a location map
close all;
cd /Users/ellynenderlin/Research/miscellaneous/RAMP
[A,S] = readgeoraster('Antarctic_RAMP_image_v2_1km.tif');
IM.x = S.XWorldLimits(1)+0.5*S.CellExtentInWorldX:S.CellExtentInWorldX:S.XWorldLimits(2)-0.5*S.CellExtentInWorldX;
IM.y = S.YWorldLimits(2)-0.5*S.CellExtentInWorldY:-S.CellExtentInWorldY:S.YWorldLimits(1)+0.5*S.CellExtentInWorldY;
IM.z=single(A);
clear A S;
%set-up a map
figureA=figure; set(gcf,'position',[450 50 800 800]);
im_cmap = colormap(gray(10001)); im_cmap(1,:) = [1 1 1];
imagesc(IM.x,IM.y,IM.z); colormap(gca,im_cmap); hold on; axis xy equal
%set-up subplots for graphs
figureB = figure; set(gcf,'position',[50 400 800 600]);
sub1b = subplot(1,2,1); sub2b = subplot(1,2,2);
figureC = figure; set(gcf,'position',[50 50 800 1000]);

%create dummy matrices to fill with meltwater flux & submerged area for E
%vs W Antarctica in order to fit trendlines to the 2 datasets
WdVdt = []; WAsub = []; EdVdt = []; EAsub = [];

%% create subplots & add site info to the map
cd(root_dir);

%set-up dummy vectors to fill with concatenated variables for all sites
start_yr = []; end_yr = []; avg_x = []; avg_y = []; depth = []; depth_uncert = []; subarea = []; subarea_uncert = [];
meltflux = []; meltflux_uncert = []; meltrate = []; meltrate_uncert = [];
%start & end coordinates for each iceberg
xcoord_o = []; ycoord_o = [];
xcoord_f = []; ycoord_f = [];

%plot
for i = 1:size(region,2)
%     cd_to_dir = ['cd ',char(region(i))]; eval(cd_to_dir);
    cd(char(region(i)));
    meltinfo = dir('*iceberg_meltinfo.csv'); landsats = dir('LC*PS.TIF');
    avgx = []; avgy = []; meltrate_v_draft = [];
%     date_o = []; xcoord_o = []; ycoord_o = [];
%     date_f = []; xcoord_f = []; ycoord_f = [];
%     flux = []; sub_area = []; meltrate = []; keeld = [];
    
    %map for study site
    figureD = figure; set(gcf,'position',[650 150 800 450]);
    [A,S] = readgeoraster(landsats(1).name);
    im.x = S.XWorldLimits(1):S.SampleSpacingInWorldX:S.XWorldLimits(2);
    im.y = S.YWorldLimits(2):-S.SampleSpacingInWorldY:S.YWorldLimits(1);
    im.z = double(A); clear A S;
    imagesc(im.x,im.y,im.z); axis xy equal; colormap(gray(10001)); hold on;
    
    %scatterplot for study site
    figureE = figure; set(gcf,'position',[1200 50 800 450]);
    
    %loop through each date pair & plot data
    for j = 1:length(meltinfo)
%         load_data = ['M=dlmread(''',meltinfo(j).name,''');']; eval(load_data);
        M=readtable(meltinfo(j).name); %table exported using plot_export_iceberg_melt_data.m
%         disp(M.Properties.VariableNames); %uncomment to display table headers prior to array conversion
        M = table2array(M); %convert to array to enable easier indexing
        
        %identify data with clear issues
        bad_data = find(M(:,18)<0); M(bad_data,:) = [];
        
        %pull variables
        dt = M(:,1); %time
        xo = M(:,2); yo = M(:,3); zo = M(:,4); rhoo = M(:,5); Vo = M(:,6); %initial locations, median elev, density, volume
        xf = M(:,7); yf = M(:,8); zf = M(:,9); rhof = M(:,10); Vf = M(:,11); %same as above but final
        coregzo = M(:,12); coregzf = M(:,13);
        dz = M(:,14); dz_sigma = M(:,15);
        dVdt = M(:,16); dVdt_uncert = M(:,17);

        %recalculate draft & submerged area to make sure they are consistent (methods may have been adjusted slightly over time)
        draft = (nanmean([rhoo rhof],2)./(repmat(rho_sw,length(nanmean([rhoo rhof],2)),1)-nanmean([rhoo rhof],2))).*nanmean([zo zf],2); %draft = M(:,18); 
        draft_uncert = M(:,19);
        Asurf = M(:,20); Asurf_uncert = M(:,21);
        lat_area = M(:,22) - Asurf; perim = lat_area./draft; clear lat_area; lat_area = perim.*draft; Asub = lat_area + Asurf; clear lat_area perim; %Asub = M(:,22); 
        Asub_uncert = M(:,23); 
        m = dVdt./Asub; %melt rate variable for plotting

%         date_o = [date_o; repmat(str2num(meltinfo(j).name(4:11)),size(xo))]; xcoord_o = [xcoord_o; xo]; ycoord_o = [ycoord_o; yo];
%         date_f = [date_f; repmat(str2num(meltinfo(j).name(13:20)),size(xf))]; xcoord_f = [xcoord_f; xf]; ycoord_f = [ycoord_f; yf];
%         flux = [flux; dVdt]; sub_area = [sub_area; Asub]; meltrate = [meltrate; (dVdt./Asub)]; keeld = [keeld; draft];

        %compile data to create a concatenated table for all sites, all dates
        yr_o = str2num(meltinfo(j).name(end-49:end-46)); 
        if mod(yr_o,4) == 0; decidate_o = (cumdays_leap(str2num(meltinfo(j).name(end-45:end-44)))+str2num(meltinfo(j).name(end-43:end-42)))/sum(modays_leap);
        else; decidate_o = (cumdays_norm(str2num(meltinfo(j).name(end-45:end-44)))+str2num(meltinfo(j).name(end-43:end-42)))/sum(modays_norm); end
        yr_f = str2num(meltinfo(j).name(end-34:end-31));
        if mod(yr_f,4) == 0; decidate_f = (cumdays_leap(str2num(meltinfo(j).name(end-30:end-29)))+str2num(meltinfo(j).name(end-28:end-27)))/sum(modays_leap);
        else; decidate_f = (cumdays_norm(str2num(meltinfo(j).name(end-30:end-29)))+str2num(meltinfo(j).name(end-28:end-27)))/sum(modays_norm); end
        start_yr = [start_yr; repmat(yr_o+decidate_o,length(draft),1)]; end_yr = [end_yr; repmat(yr_f+decidate_f,length(draft),1)];
        plot_yrs = [yr_o+decidate_o yr_f+decidate_f];
        clear yr_* decidate_*;
        avgx = [avgx; nanmean([xo xf],2)]; avgy = [avgy; nanmean([yo yf],2)]; %average coordinates for regional map
        xcoord_o = [xcoord_o; xo]; ycoord_o = [ycoord_o; yo]; xcoord_f = [xcoord_f; xf]; ycoord_f = [ycoord_f; yf]; %coordinates for coordinate data table
        avg_x = [avg_x; nanmean([xo xf],2)]; avg_y = [avg_y; nanmean([yo yf],2)]; %average coordinates for site map & regional data table
        depth = [depth; draft]; depth_uncert = [depth_uncert; draft_uncert]; %keel depth
        subarea = [subarea; Asub]; subarea_uncert = [subarea_uncert; Asub_uncert]; %submerged area
        meltflux = [meltflux; dVdt]; meltflux_uncert = [meltflux_uncert; dVdt_uncert]; %meltwater flux
        meltrate = [meltrate; m]; meltrate_uncert = [meltrate_uncert; abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2)]; %melt rate
        
        
        %multi-panel subplots of all data
        figure(figureB);
        subplot(sub1b);
        errorbar(Asub,dVdt/86400,dVdt_uncert/86400,dVdt_uncert/86400,Asub_uncert,Asub_uncert,plot_marker,...
            'markerfacecolor','w','markersize',symbol_size,'markeredgecolor',region_colors(i,:),...
            'color',region_colors(i,:)); hold on;
        pl(i) = plot(Asub,dVdt/86400,plot_marker,'markerfacecolor',region_colors(i,:),'markerfacecolor',region_colors(i,:),...
            'color',region_colors(i,:),'markersize',symbol_size,'linewidth',1); hold on;
        subplot(sub2b);
        errorbar(draft,365*m,365*abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),365*abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,plot_marker,...
            'markerfacecolor','w','markersize',symbol_size,'markeredgecolor',region_colors(i,:),...
            'color',region_colors(i,:)); hold on;
        plot(draft,365*m,plot_marker,'markerfacecolor',region_colors(i,:),'markerfacecolor',region_colors(i,:),...
            'color',region_colors(i,:),'markersize',symbol_size,'linewidth',1); hold on;
        
        %subplots for each region but all in one figure window
        figure(figureC);
        if j == 1
            subpl = subplot(8,2,plot_loc(i));
        else
            subplot(subpl);
        end
        %add dummy data for the legend
        if plot_loc(i) == 14
            if j == 1
                for k = 1:length(yrs)
                    yrpl(k) = plot(draft(1),365*dVdt(1)./Asub(1),plot_marker,'markerfacecolor',plot_color(yrs(k)-(min(yrs)-1),:),...
                        'markeredgecolor',plot_color(yrs(k)-(min(yrs)-1),:),'markersize',symbol_size,...
                        'color',plot_color(yrs(k)-(min(yrs)-1),:)); hold on;
                end
            end
        end
        %adjust subplot positions
        plotpos = get(gca,'position');
        if j == 1
            if mod(plot_loc(i),2) ~= 0
                set(gca,'position',[plotpos(1)-0.03 plotpos(2)+0.04/(plot_loc(i)+1) 1.1*plotpos(3) 1.1*plotpos(4)]);
            else
                set(gca,'position',[plotpos(1) plotpos(2)+0.04/(plot_loc(i)) 1.1*plotpos(3) 1.1*plotpos(4)]);
            end
        end
        clear plotpos;
        %plot face color color-coded by date
        errorbar(draft,365*m,365*abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),365*abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,[plot_marker,'k'],...
            'markerfacecolor','none','markeredgecolor',plot_color(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markersize',symbol_size,...
            'color',plot_color(round(nanmean(plot_yrs))-(min(yrs)-1),:)); hold on;
%         scatter(draft,365*m,6*symbol_size,plot_color(round(nanmean(plot_yrs))-(min(yrs)-1),:),plot_marker); hold on;
        plot(draft,365*m,plot_marker,'markerfacecolor',plot_color(round(nanmean(plot_yrs))-(min(yrs)-1),:),...
            'markerfacecolor',plot_color(round(nanmean(plot_yrs))-(min(yrs)-1),:),...
            'color',plot_color(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markersize',symbol_size,'linewidth',1); hold on;
        
        %add data to the site map-view figure
        figure(figureD);
        for k = 1:length(m)
            symbol_color = ceil(2*(m(k))*365);
            if symbol_color > length(highmelt_cmap); symbol_color = length(highmelt_cmap); end
            plot(nanmean([xo(k) xf(k)]),nanmean([yo(k) yf(k)]),[map_marker,'k'],'markerfacecolor',highmelt_cmap(symbol_color,:),'markersize',round(draft(k)/30)+10); hold on;
            clear symbol_color;
        end
        
        %add data to site scatterplots
        figure(figureE);
        errorbar(draft,365*m,365*abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),365*abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,[plot_marker,'k'],...
            'markerfacecolor','none','markeredgecolor',plot_color(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markersize',2*symbol_size,...
            'color',plot_color(round(nanmean(plot_yrs))-(min(yrs)-1),:)); hold on;
%         indpl(j) = scatter(draft,365*m,2*(6*symbol_size),plot_color(round(nanmean(plot_yrs))-(min(yrs)-1),:),plot_marker); hold on;
        indpl(j) = plot(draft,365*m,plot_marker,'markerfacecolor',plot_color(round(nanmean(plot_yrs))-(min(yrs)-1),:),...
            'markerfacecolor',plot_color(round(nanmean(plot_yrs))-(min(yrs)-1),:),...
            'color',plot_color(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markersize',2*symbol_size,'linewidth',1); hold on;
        daterange(j,:) = [meltinfo(j).name(end-49:end-46),'/',meltinfo(j).name(end-45:end-44),'/',meltinfo(j).name(end-43:end-42),'-',meltinfo(j).name(end-34:end-31),'/',meltinfo(j).name(end-30:end-29),'/',meltinfo(j).name(end-28:end-27)];
        
        %remove date-specific variables
        clear m dt xo yo zo rhoo Vo xf yf zf rhof Vf coregzo coregzf dz dz_sigma dVdt dVdt_uncert draft draft_uncert Asurf Asurf_uncert Asub Asub_uncert;
    end
    
    %add a symbol to the site map
    figure(figureA); colormap(gca,'gray');
    mp(i) = plot(nanmean(avgx),nanmean(avgy),[map_marker,'k'],'markerfacecolor',region_colors(i,:),'markeredgecolor','k','markersize',12); hold on;
    
    %format the scatterplot subplots
    figure(figureC);
    if nanmean(avgx) > -1.5e6
        set(gca,'xlim',[0 400],'xtick',[0:100:400],'xticklabel',[],...
            'ylim',[0 9],'ytick',[0:4:9],'yticklabel',[0:4:9],'fontsize',16); grid on;
    else
        if nanmean(avgx) > -2.45e6 && nanmean(avgy) > 1.2e6
            set(gca,'xlim',[0 400],'xtick',[0:100:400],'xticklabel',[],...
                'ylim',[0 22],'ytick',[0:10:22],'yticklabel',[0:10:22],'fontsize',16); grid on;
        else
            set(gca,'xlim',[0 800],'xtick',[0:200:800],'xticklabel',[],...
                'ylim',[0 70],'ytick',[0:30:70],'yticklabel',[0:30:70],'fontsize',16); grid on;
        end
    end
    %add axes labels
    if plot_loc(i) == 14
        legsub = legend(yrpl,num2str(yrs')); set(legsub,'location','southoutside','NumColumns',5,'fontsize',16);
        legpos = get(legsub,'position'); set(legsub,'position',[0.5 0.13 legpos(3) legpos(4)]); clear legpos;
        set(gca,'xlim',[0 400],'xtick',[0:100:400],'xticklabel',[0:100:400]); grid on;
        xlabel('Draft (m b.s.l.)','fontsize',16);
    end
    if plot_loc(i) == 15
        xlabel('Draft (m b.s.l.)','fontsize',16); ylb = ylabel('Iceberg melt rate (m yr^{-1})','fontsize',16);
        set(gca,'xlim',[0 800],'xtick',[0:200:800],'xticklabel',[0:200:800]);
        set(ylb,'position',[-100 350 -1]);
    end
    text(0.6*max(get(gca,'xlim')),0.85*max(get(gca,'ylim')),char(disp_names(i)),'fontsize',16);
    
    %format the site map-view figure
    figure(figureD);
    if sqrt((max(avgx)-min(avgx)).^2 + (max(avgy)-min(avgy)).^2)+10000 < 50000
        set(gca,'xlim',[min(avgx)-5000 max(avgx)+5000],'xtick',[(ceil(min(avgx)/1000)*1000-5000):5000:(floor(max(avgx)/1000)*1000+5000)],'xticklabel',[(ceil(min(avgx)/1000)-5):5:(floor(max(avgx)/1000)+5)],...
            'ylim',[min(avgy)-5000 max(avgy)+5000],'ytick',[(ceil(min(avgy)/1000)*1000-5000):5000:(floor(max(avgy)/1000)*1000+5000)],'yticklabel',[(ceil(min(avgy)/1000)-5):5:(floor(max(avgy)/1000)+5)],...
            'fontsize',16);
    else
        set(gca,'xlim',[min(avgx)-5000 max(avgx)+5000],'xtick',[(ceil(min(avgx)/1000)*1000-5000):10000:(floor(max(avgx)/1000)*1000+5000)],'xticklabel',[(ceil(min(avgx)/1000)-5):10:(floor(max(avgx)/1000)+5)],...
            'ylim',[min(avgy)-5000 max(avgy)+5000],'ytick',[(ceil(min(avgy)/1000)*1000-5000):5000:(floor(max(avgy)/1000)*1000+5000)],'yticklabel',[(ceil(min(avgy)/1000)-5):5:(floor(max(avgy)/1000)+5)],...
            'fontsize',16);
    end
    xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
    %resize vertical dimension to maximize figure window usage
    xlims = get(gca,'xlim'); ylims = get(gca,'ylim'); figpos = get(gcf,'position');
    set(gcf,'position',[figpos(1) figpos(2) figpos(3) (max(ylims)-min(ylims))/(max(xlims)-min(xlims))*figpos(3)]);
    %add legends
    if ~isempty(strmatch('Thwaites',char(region(i))))
        rectangle('position',[min(xlims)+0.80*(max(xlims)-min(xlims)) min(ylims)+0.725*(max(ylims)-min(ylims)) 0.15*(max(xlims)-min(xlims)) 0.175*(max(ylims)-min(ylims))],'curvature',[0,0],'facecolor','w','linewidth',2);
        plot(min(xlims)+0.825*(max(xlims)-min(xlims)),min(ylims)+0.87*(max(ylims)-min(ylims)),[map_marker,'k'],'markersize',round(50/30+10),'markerfacecolor','w'); text(min(xlims)+0.85*(max(xlims)-min(xlims)),min(ylims)+0.8675*(max(ylims)-min(ylims)),'50 m','fontsize',16);
        plot(min(xlims)+0.825*(max(xlims)-min(xlims)),min(ylims)+0.82*(max(ylims)-min(ylims)),[map_marker,'k'],'markersize',round(150/30+10),'markerfacecolor','w'); text(min(xlims)+0.85*(max(xlims)-min(xlims)),min(ylims)+0.8175*(max(ylims)-min(ylims)),'150 m','fontsize',16);
        plot(min(xlims)+0.825*(max(xlims)-min(xlims)),min(ylims)+0.76*(max(ylims)-min(ylims)),[map_marker,'k'],'markersize',round(300/30+10),'markerfacecolor','w'); text(min(xlims)+0.85*(max(xlims)-min(xlims)),min(ylims)+0.7575*(max(ylims)-min(ylims)),'300 m','fontsize',16);
        rectangle('position',[min(xlims)+0.575*(max(xlims)-min(xlims)) min(ylims)+0.725*(max(ylims)-min(ylims)) 0.20*(max(xlims)-min(xlims)) 0.25*(max(ylims)-min(ylims))],'curvature',[0,0],'facecolor','w','linewidth',2);
        %             symbol_color = round(([0.001 0.01 0.1])*1000)+1; symbol_color(symbol_color>length(highmelt_cmap)) = length(highmelt_cmap);
        for k = 1:length(highmelt_cmap)
            plot([min(xlims)+0.60*(max(xlims)-min(xlims)) min(xlims)+0.64*(max(xlims)-min(xlims))],...
                [min(ylims)+0.945*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap)) min(ylims)+0.945*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap))],'-','linewidth',2*((max(ylims)-min(ylims))/(max(xlims)-min(xlims))),'color',highmelt_cmap(k,:));
        end
        text(min(xlims)+0.65*(max(xlims)-min(xlims)),min(ylims)+0.95*(max(ylims)-min(ylims))-((0.004*(max(ylims)-min(ylims)))),'1 m yr^{-1}','fontsize',16);
        text(min(xlims)+0.65*(max(xlims)-min(xlims)),min(ylims)+0.95*(max(ylims)-min(ylims))-((0.04*(max(ylims)-min(ylims)))),'10 m yr^{-1}','fontsize',16);
        text(min(xlims)+0.65*(max(xlims)-min(xlims)),min(ylims)+0.95*(max(ylims)-min(ylims))-((0.20*(max(ylims)-min(ylims)))),'50 m yr^{-1}','fontsize',16);
    elseif ~isempty(strmatch('Mertz',char(region(i))))
        rectangle('position',[min(xlims)+0.25*(max(xlims)-min(xlims)) min(ylims)+0.725*(max(ylims)-min(ylims)) 0.15*(max(xlims)-min(xlims)) 0.175*(max(ylims)-min(ylims))],'curvature',[0,0],'facecolor','w','linewidth',2);
        plot(min(xlims)+0.275*(max(xlims)-min(xlims)),min(ylims)+0.87*(max(ylims)-min(ylims)),[map_marker,'k'],'markersize',round(50/30+10),'markerfacecolor','w'); text(min(xlims)+0.30*(max(xlims)-min(xlims)),min(ylims)+0.8675*(max(ylims)-min(ylims)),'50 m','fontsize',16);
        plot(min(xlims)+0.275*(max(xlims)-min(xlims)),min(ylims)+0.82*(max(ylims)-min(ylims)),[map_marker,'k'],'markersize',round(150/30+10),'markerfacecolor','w'); text(min(xlims)+0.30*(max(xlims)-min(xlims)),min(ylims)+0.8175*(max(ylims)-min(ylims)),'150 m','fontsize',16);
        plot(min(xlims)+0.275*(max(xlims)-min(xlims)),min(ylims)+0.76*(max(ylims)-min(ylims)),[map_marker,'k'],'markersize',round(300/30+10),'markerfacecolor','w'); text(min(xlims)+0.30*(max(xlims)-min(xlims)),min(ylims)+0.7575*(max(ylims)-min(ylims)),'300 m','fontsize',16);
        rectangle('position',[min(xlims)+0.025*(max(xlims)-min(xlims)) min(ylims)+0.725*(max(ylims)-min(ylims)) 0.20*(max(xlims)-min(xlims)) 0.25*(max(ylims)-min(ylims))],'curvature',[0,0],'facecolor','w','linewidth',2);
        %             symbol_color = round(([0.001 0.01 0.1])*1000)+1; symbol_color(symbol_color>length(highmelt_cmap)) = length(highmelt_cmap);
        for k = 1:length(highmelt_cmap)
            plot([min(xlims)+0.05*(max(xlims)-min(xlims)) min(xlims)+0.09*(max(xlims)-min(xlims))],...
                [min(ylims)+0.945*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap)) min(ylims)+0.945*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap))],'-','linewidth',2*((max(ylims)-min(ylims))/(max(xlims)-min(xlims))),'color',highmelt_cmap(k,:));
        end
        text(min(xlims)+0.10*(max(xlims)-min(xlims)),min(ylims)+0.95*(max(ylims)-min(ylims))-((0.004*(max(ylims)-min(ylims)))),'1 m yr^{-1}','fontsize',16);
        text(min(xlims)+0.10*(max(xlims)-min(xlims)),min(ylims)+0.95*(max(ylims)-min(ylims))-((0.04*(max(ylims)-min(ylims)))),'10 m yr^{-1}','fontsize',16);
        text(min(xlims)+0.10*(max(xlims)-min(xlims)),min(ylims)+0.95*(max(ylims)-min(ylims))-((0.20*(max(ylims)-min(ylims)))),'50 m yr^{-1}','fontsize',16);
    else
        rectangle('position',[min(xlims)+0.25*(max(xlims)-min(xlims)) min(ylims)+0.020*(max(ylims)-min(ylims)) 0.15*(max(xlims)-min(xlims)) 2800],'curvature',[0,0],'facecolor','w','linewidth',2); %scaling height = 0.175*(max(ylims)-min(ylims))
        plot(min(xlims)+0.275*(max(xlims)-min(xlims)),min(ylims)+0.020*(max(ylims)-min(ylims))+2260,[map_marker,'k'],'markersize',round(50/30+10),'markerfacecolor','w'); text(min(xlims)+0.30*(max(xlims)-min(xlims)),min(ylims)+0.020*(max(ylims)-min(ylims))+2260,'50 m','fontsize',16); %scaling y-offset = 0.16*(max(ylims)-min(ylims))
        plot(min(xlims)+0.275*(max(xlims)-min(xlims)),min(ylims)+0.020*(max(ylims)-min(ylims))+1460,[map_marker,'k'],'markersize',round(150/30+10),'markerfacecolor','w'); text(min(xlims)+0.30*(max(xlims)-min(xlims)),min(ylims)+0.020*(max(ylims)-min(ylims))+1460,'150 m','fontsize',16); %scaling y-offset = 0.11*(max(ylims)-min(ylims))
        plot(min(xlims)+0.275*(max(xlims)-min(xlims)),min(ylims)+0.020*(max(ylims)-min(ylims))+500,[map_marker,'k'],'markersize',round(300/30+10),'markerfacecolor','w'); text(min(xlims)+0.30*(max(xlims)-min(xlims)),min(ylims)+0.020*(max(ylims)-min(ylims))+500,'300 m','fontsize',16); %scaling y-offset = 0.05*(max(ylims)-min(ylims))
        rectangle('position',[min(xlims)+0.025*(max(xlims)-min(xlims)) min(ylims)+0.020*(max(ylims)-min(ylims)) 0.20*(max(xlims)-min(xlims)) 0.255*(max(ylims)-min(ylims))],'curvature',[0,0],'facecolor','w','linewidth',2);
        %             symbol_color = round(([0.001 0.01 0.1])*1000)+1; symbol_color(symbol_color>length(highmelt_cmap)) = length(highmelt_cmap);
        for k = 1:length(highmelt_cmap)
            plot([min(xlims)+0.05*(max(xlims)-min(xlims)) min(xlims)+0.09*(max(xlims)-min(xlims))],...
                [min(ylims)+0.245*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap)) min(ylims)+0.245*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap))],'-','linewidth',2*((max(ylims)-min(ylims))/(max(xlims)-min(xlims))),'color',highmelt_cmap(k,:));
        end
        text(min(xlims)+0.10*(max(xlims)-min(xlims)),min(ylims)+0.25*(max(ylims)-min(ylims))-((0.004*(max(ylims)-min(ylims)))),'1 m yr^{-1}','fontsize',16);
        text(min(xlims)+0.10*(max(xlims)-min(xlims)),min(ylims)+0.25*(max(ylims)-min(ylims))-((0.04*(max(ylims)-min(ylims)))),'10 m yr^{-1}','fontsize',16);
        text(min(xlims)+0.10*(max(xlims)-min(xlims)),min(ylims)+0.25*(max(ylims)-min(ylims))-((0.20*(max(ylims)-min(ylims)))),'50 m yr^{-1}','fontsize',16);
    end
    clear im;
    saveas(figureD,[char(region(i)),'_melt-map.eps'],'epsc'); saveas(figureD,[char(region(i)),'_melt-map.png'],'png');
    clear xlims ylims;
    
    %format the site scatterplot
    figure(figureE);
    xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
    if max(xlims) > 500
        set(gca,'xlim',[0 max(xlims)],'xtick',[0:100:round(max(xlims))],'xticklabel',[0:100:round(max(xlims))],...
            'ylim',[0 max(ylims)],'ytick',[0:round(max(ylims)/3):max(ylims)],'yticklabel',[0:round(max(ylims)/3):max(ylims)],'fontsize',16); grid on;
    else
        set(gca,'xlim',[0 max(xlims)],'xtick',[0:50:round(max(xlims))],'xticklabel',[0:50:round(max(xlims))],...
            'ylim',[0 max(ylims)],'ytick',[0:round(max(ylims)/3):max(ylims)],'yticklabel',[0:round(max(ylims)/3):max(ylims)],'fontsize',16); grid on;
    end
    xlabel('Draft (m b.s.l.)','fontsize',16); ylabel('Melt rate (m yr^{-1})','fontsize',16);
    %find the best location for the legend based on data
    leg = legend(indpl,daterange); 
    if i <= 2
        set(leg,'fontsize',16,'location','southeast');
    else
        set(leg,'fontsize',16,'location','northwest');
    end
    saveas(figureE,[char(region(i)),'_iceberg-meltrate-v-draft.eps'],'epsc'); saveas(figureE,[char(region(i)),'_iceberg-meltrate-v-draft.png'],'png');
    clear xlims ylims;
    
    %save all dates and coords to a tab-delimited text file
    cd ..
    clear subpl indpl daterange;
end
figure(figureA);
set(gca,'xlim',[-28e5 28e5],'xtick',[-32e5:8e5:32e5],'xticklabel',[-3200:800:3200],...
    'ylim',[-24e5 24e5],'ytick',[-24e5:8e5:24e5],'yticklabel',[-2400:800:2400],'fontsize',16); grid off;
xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
graticuleps(-50:-5:-90,-180:30:180);
text(0,6.5e5,['85',char(176),'S'],'fontsize',16); text(0,12.0e5,['80',char(176),'S'],'fontsize',16); 
text(0,17.5e5,['75',char(176),'S'],'fontsize',16); text(0,23.0e5,['70',char(176),'S'],'fontsize',16);
text(-16.5e5,25.25e5,['-30',char(176),'E'],'fontsize',16); text(12.5e5,25.25e5,['30',char(176),'E'],'fontsize',16); 
colormap(gca,im_cmap);
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
legmap = legend(mp,[char(leg_names)]); set(legmap,'location','westoutside','fontsize',16); 
saveas(gcf,'Antarctic-iceberg-map.eps','epsc'); saveas(gcf,'Antarctic-iceberg-map.png','png');

%save the subplots containing all data
figure(figureB);
subplot(sub1b);
leg1 = legend(pl,[char(leg_names)]); set(leg1,'location','westoutside','fontsize',16); %set(leg1,'position',[0.03    0.2756    0.0685    0.4838]);
set(gca,'xlim',[0 7e6],'xtick',[0:1e6:7e6],'xticklabel',[0:1:7],...
    'ylim',[0 6.8],'ytick',[0:1:7],'yticklabel',[0:1:7],'fontsize',16); grid on;
xlabel('Submerged area (km^2)','fontsize',16); ylabel('Meltwater flux (m^3 s^{-1})','fontsize',16);
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
text(0.05*max(xlims),0.95*max(ylims),'a) ','color','k','fontsize',16);
subplot(sub2b);
set(gca,'xlim',[0 800],'xtick',[0:200:800],'xticklabel',[0:200:800],...
    'ylim',[0 72],'ytick',[0:24:72],'yticklabel',[0:24:72],'fontsize',16); grid on;
xlabel('Draft (m b.s.l.)','fontsize',16); ylabel('Melt rate (m yr^{-1})','fontsize',16);
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
text(0.05*max(xlims),0.95*max(ylims),'b) ','color','k','fontsize',16);
saveas(gcf,'Antarctic-iceberg-lumped-plots.eps','epsc'); saveas(gcf,'Antarctic-iceberg-lumped-plots.png','png');

%save the subplots sorted by study site
figure(figureC);
saveas(gcf,'Antarctic-iceberg-subplots.eps','epsc'); saveas(gcf,'Antarctic-iceberg-subplots.png','png');


%save the data tables
%MELT
column_names = {'Start Date' 'End Date' 'Polar Stereo Easting' 'Polar Stereo Northing'...
    'Average Median Draft' 'Median Draft Variability' 'Average Submerged Area' 'Submerged Area Variability'...
    'Meltwater Flux' 'Meltwater Flux Uncertainty' 'Melt Rate' 'Melt Rate Uncertainty'};
column_units = {'years' 'years' 'meters' 'meters'...
    'meters b.s.l.' 'meters b.s.l.' 'cubic meters' 'cubic meters'...
    'cubic meters per day' 'cubic meters per day' 'meters per day' 'meters per day'};
T=table(start_yr,end_yr,avg_x,avg_y,depth,depth_uncert,subarea,subarea_uncert,meltflux,meltflux_uncert,meltrate,meltrate_uncert);
T.Properties.VariableNames = column_names; T.Properties.VariableUnits = column_units;
writetable(T,[root_dir,'Antarctic-iceberg-meltinfo.csv']);
disp('Antarctic iceberg melt rate text file written');
clear column_*;
%COORDINATES ONLY
column_names = {'Start Date' 'Start Polar Stereo Easting' 'Start Polar Stereo Northing'...
    'End Date' 'End Polar Stereo Easting' 'End Polar Stereo Northing'};
column_units = {'years' 'meters' 'meters'...
    'years' 'meters' 'meters'};
t=table(start_yr,xcoord_o,ycoord_o,end_yr,xcoord_f,ycoord_f);
t.Properties.VariableNames = column_names; t.Properties.VariableUnits = column_units;
writetable(t,[root_dir,'Antarctic-iceberg-PScoords.csv']);
disp('Antarctic iceberg coordinates text file written');
clear T t;

%% replot Larsen A & B to hone-in on potential variations btw dates
% cd /users/mariamadryak/Desktop/Antarctic_icebergs
% clear marker; yrs = [2011:1:2018]; 
% region = [{'LarsenA'},{'LarsenB'},{'Bugge'},{'Marguerite'},{'Biscoe'},{'Palmer'},{'Danco'}];
% marker = ['x','*','s','d','^','p','h','o'];
% plot_cmap = colormap(parula(length(yrs))); %colormap for emphasizing different locations
% for i = 1:2
%     figureB = figure; set(gcf,'position',[50 400 2000 800]);
%     cd_to_dir = ['cd ',char(region(i))]; eval(cd_to_dir);
%     meltinfo = dir('*iceberg_meltinfo.txt'); landsats = dir('LC*PS.TIF');
%     xcoord_o = []; ycoord_o = []; flux = []; sub_area = []; meltrate = []; keeld = [];
%     for j = 1:length(meltinfo)
%         plot_yrs = [str2num(meltinfo(j).name(end-37:end-34)) str2num(meltinfo(j).name(end-28:end-25))];
%         load_data = ['M=dlmread(''',meltinfo(j).name,''');']; eval(load_data);
%         dt = M(:,1); %time
%         xo = M(:,2); yo = M(:,3); zo = M(:,4); rhoo = M(:,5); Vo = M(:,6); %initial locations, median elev, density, volume
%         xf = M(:,7); yf = M(:,8); zf = M(:,9); rhof = M(:,10); Vf = M(:,11); %same as above but final
%         coregzo = M(:,12); coregzf = M(:,13);
%         dz = M(:,14); dz_sigma = M(:,15);
%         dVdt = M(:,16); dVdt_uncert = M(:,17);
%         draft = M(:,18); draft_uncert = M(:,19);
%         Asurf = M(:,20); Asurf_uncert = M(:,21);
%         Asub = M(:,22); Asub_uncert = M(:,23);
%         avgx(i) = nanmean(xo); avgy(i) = nanmean(yo); %this determines where symbols get plotted (i.e. study site locations)
%         %     region(i,:) = 'Larsen A  ';
%         xcoord_o = [xcoord_o xo']; ycoord_o = [ycoord_o yo'];
%         flux = [flux dVdt']; sub_area = [sub_area Asub']; meltrate = [meltrate (dVdt./Asub)']; keeld = [keeld draft'];
%         if ~isempty(strmatch(marker(i),'p')) || ~isempty(strmatch(marker(i),'h')); symbol_size = 14; else symbol_size = 12; end
%         
%         figure(figureB);
%         errorbar(draft,dVdt./Asub,abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,[marker(round(nanmean(plot_yrs))-(min(yrs)-1)),'k'],'markerfacecolor',plot_cmap(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markersize',symbol_size); hold on;
%         plot(draft,dVdt./Asub,[marker(round(nanmean(plot_yrs))-(min(yrs)-1)),'k'],'markerfacecolor',plot_cmap(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markersize',symbol_size); hold on;
%         
%         clear dt xo yo zo rhoo Vo xf yf zf rhof Vf coregzo coregzf dz dz_sigma dVdt dVdt_uncert draft draft_uncert Asurf Asurf_uncert Asub Asub_uncert;
%     end
%     set(gca,'xlim',[0 300],'xtick',[0:100:300],'xticklabel',[0:100:300],...
%                 'ylim',[0 0.04],'ytick',[0:0.01:0.04],'yticklabel',[0:1:4],'fontsize',16); grid on;
%     text(0.1*max(get(gca,'xlim')),0.9*max(get(gca,'ylim')),char(region(i)),'fontsize',16);
% %     %add a map-view figure
% %     if ~isempty(landsats)
% %         [A,S] = geotiffread(landsats(1).name);
% %         im.x = S.XWorldLimits(1):S.SampleSpacingInWorldX:S.XWorldLimits(2);
% %         im.y = S.YWorldLimits(2):-S.SampleSpacingInWorldY:S.YWorldLimits(1);
% %         im.z = double(A); clear A S;
% %         figureD = figure; set(gcf,'position',[450 450 900 500]);
% %         imagesc(im.x,im.y,im.z); axis xy equal; colormap gray; hold on;
% %         for j = 1:length(meltrate)
% %             symbol_color = round((meltrate(j)-min(meltrate))*1000)+1;
% %             if i>=8
% %                 if symbol_color > length(highmelt_cmap); symbol_color = length(highmelt_cmap); end
% %                 plot(xcoord_o(j),ycoord_o(j),'ok','markerfacecolor',highmelt_cmap(symbol_color,:),'markersize',round(keeld(j)/30)+10); hold on;
% %             else
% %                 if symbol_color > length(lowmelt_cmap); symbol_color = length(lowmelt_cmap); end
% %                 plot(xcoord_o(j),ycoord_o(j),'ok','markerfacecolor',lowmelt_cmap(symbol_color,:),'markersize',round(keeld(j)/30)+10); hold on;
% %             end
% %             clear symbol_color;
% %         end
% %         set(gca,'xlim',[min(xcoord_o)-5000 max(xcoord_o)+5000],'xtick',[(ceil(min(xcoord_o)/1000)*1000-5000):5000:(floor(max(xcoord_o)/1000)*1000+5000)],'xticklabel',[(ceil(min(xcoord_o)/1000)-5):5:(floor(max(xcoord_o)/1000)+5)],...
% %             'ylim',[min(ycoord_o)-5000 max(ycoord_o)+5000],'ytick',[(ceil(min(ycoord_o)/1000)*1000-5000):5000:(floor(max(ycoord_o)/1000)*1000+5000)],'yticklabel',[(ceil(min(ycoord_o)/1000)-5):5:(floor(max(ycoord_o)/1000)+5)],...
% %             'fontsize',16);
% %         xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
% %         clear im;
% %         saveas(figureD,[char(region(i)),'_melt-map.eps'],'epsc'); saveas(figureD,[char(region(i)),'_melt-map.png'],'png');
% %     end
%     cd ..
% end