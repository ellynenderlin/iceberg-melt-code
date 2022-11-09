%PLOT_MELTRATE_FIGS_UPDATED: Generate maps and melt rate subplots for all
%Antarctic iceberg melt datasets. Includes ocean observation data in
%figures. Exports concatenated data for all sites to a table.


%% Initialize
clearvars; close all; drawnow;
addpath('/users/ellynenderlin/Research/miscellaneous/general-code','/users/ellynenderlin/Research/miscellaneous/general-code/cmocean');

%specify paths & file names for data
iceberg_path = '/Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt/';
figure_path = [iceberg_path,'figures/'];

%specify generic variables
plot_yrs = []; avgx = []; avgy = []; region = []; warning off;
rho_sw = 1026; %sea water density in kg m^-3
years = [2011.75 2022.25]; year_ticks = [2013:2:2022]; %approximate date range for plots

%specify plot parameters
yrs = [round(min(years)):1:round(max(years))]; 
plot_color = cmocean('tempo',length(yrs)+1); plot_color = plot_color(2:end,:); %colormap for emphasizing different years
highmelt_cmap = cmocean('amp',100); % highmelt_cmap = flipud(colormap(hot(100)));
% plot_cmap = colormap(parula(15)); %colormap for emphasizing different locations
region = [{'Edgeworth-LarsenA'},{'Crane-LarsenB'},{'Ronne'},{'Filchner'},{'Amery'},{'Totten'},{'Mertz'},...
    {'Thwaites'},{'Ferrigno-Eltanin'},{'Seller-Bugge'},{'Heim-Marguerite'},{'Widdowson-Biscoe'},{'Cadman-Palmer'},{'Blanchard-Danco'},{'Leonardo-Danco'}]; %site directory names
plot_letters = [{'i)'},{'j)'},{'k)'},{'l)'},{'m)'},{'n)'},{'o)'},{'h)'},{'g)'},{'f)'},{'e)'},{'d)'},{'c)'},{'b)'},{'a)'}]; %plot letters for sites to be used in geographically-arranged subplots
leg_names = [{'Edgeworth'},{'Crane'},{'Ronne'},{'Filchner'},{'Polar Times'},{'Totten'},{'Mertz'},...
    {'Thwaites'},{'Ferrigno'},{'Seller'},{'Heim'},{'Widdowson'},{'Cadman'},{'Blanchard'},{'Leonardo'}]; %generic figure labels
leg_abbrevs = [{'LA'},{'LB'},{'RI'},{'FI'},{'PT'},{'TI'},{'MI'},...
    {'TG'},{'FG'},{'SG'},{'HG'},{'WG'},{'CG'},{'BG'},{'LG'}]; %generic figure labels
leg_ref = [9,10,11,12,13,14,15,8,7,6,5,4,3,2,1]; %arrange legend in alphabetical order
region_ref = [1,1,2,2,2,2,2,3,3,4,4,4,4,4,4]; %group according to region (1=EAP, 2=EAIS, 3=WAIS, 4=WAP)

%specify plot parameters
map_marker = 's'; %marker symbol for maps
plot_marker = 's'; %marker symbol for scatterplots
symbol_size = 4; %symbols size
plot_loc = [2,4,6,8,10,12,14,15,13,11,9,7,5,3,1]; %specifies the subplot locations for the study sites listed as "region", "disp_names", and "leg_names" above
EAP_color = [77,172,38]./255; EAIS_color = [184,225,134]./255; WAIS_color = [241,182,218]./255; WAP_color = [208,28,139]./255; 
%create color pallettes
region_colors = zeros(length(region_ref),3);
region_colors(find(region_ref == 1),:) = repmat(EAP_color,length(find(region_ref == 1)),1); 
region_colors(find(region_ref == 2),:) = repmat(EAIS_color,length(find(region_ref == 2)),1); 
region_colors(find(region_ref == 3),:) = repmat(WAIS_color,length(find(region_ref == 3)),1); 
region_colors(find(region_ref == 4),:) = repmat(WAP_color,length(find(region_ref == 4)),1); 

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

close all;

%% set-up plots
close all; drawnow;

%set-up location maps: (1) alphabetical labeling & (2) site abbreviation labeling
%(1) alphabetical labeling
figureA1=figure; set(gcf,'position',[450 50 800 800]);
im_cmap = colormap(gray(10001)); im_cmap(1,:) = [1 1 1];
imagesc(IM.x,IM.y,IM.z); colormap(gca,im_cmap); hold on; axis xy equal
xlims = [-26.5e5 27.5e5]; ylims = [-22.5e5 22.5e5];
set(gca,'xlim',xlims,'xtick',[-30e5:6e5:30e5],'xticklabel',[-3000:600:3000],...
    'ylim',ylims,'ytick',[-24e5:6e5:24e5],'yticklabel',[-2400:600:2400],'fontsize',16); grid off;
xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
%add polar stereo coordinates
graticuleps(-50:-5:-90,-180:30:180);
text(-5.5e5,3.25e5,['85',char(176),'S'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',45); 
text(-9.5e5,7.0e5,['80',char(176),'S'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',45); 
text(-13.5e5,10.75e5,['75',char(176),'S'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',45); 
text(-17.55e5,14.55e5,['70',char(176),'S'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',45);
text(-21.65e5,18.4e5,['65',char(176),'S'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',45);
text(0.75e5,max(ylims)-2.75e5,['0',char(176),'E'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',90); 
text(max(xlims)-4e5,0.75e5,['90',char(176),'E'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',0);
text(0.75e5,min(ylims)+4.25e5,['180',char(176),'E'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',-90);
text(min(xlims)+0.5e5,0.75e5,['-90',char(176),'E'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',0);
text(-0.75e5,min(ylims),['-180',char(176),'E'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',90);
clear xlims ylims;
%(2) site abbreviation labelling (for USAP-DC overview fig)
figureA2=figure; set(gcf,'position',[450 50 800 800]);
im_cmap = colormap(gray(10001)); im_cmap(1,:) = [1 1 1];
imagesc(IM.x,IM.y,IM.z); colormap(gca,im_cmap); hold on; axis xy equal
xlims = [-26.5e5 27.5e5]; ylims = [-22.5e5 22.5e5];
set(gca,'xlim',xlims,'xtick',[-30e5:6e5:30e5],'xticklabel',[-3000:600:3000],...
    'ylim',ylims,'ytick',[-24e5:6e5:24e5],'yticklabel',[-2400:600:2400],'fontsize',16); grid off;
xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
%add polar stereo coordinates
graticuleps(-50:-5:-90,-180:30:180);
text(-5.5e5,3.25e5,['85',char(176),'S'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',45); 
text(-9.5e5,7.0e5,['80',char(176),'S'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',45); 
text(-13.5e5,10.75e5,['75',char(176),'S'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',45); 
text(-17.55e5,14.55e5,['70',char(176),'S'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',45);
text(-21.65e5,18.4e5,['65',char(176),'S'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',45);
text(0.75e5,max(ylims)-2.75e5,['0',char(176),'E'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',90); 
text(max(xlims)-4e5,0.75e5,['90',char(176),'E'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',0);
text(0.75e5,min(ylims)+4.25e5,['180',char(176),'E'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',-90);
text(min(xlims)+0.5e5,0.75e5,['-90',char(176),'E'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',0);
text(-0.75e5,min(ylims),['-180',char(176),'E'],'fontsize',16,'color',[0.5 0.5 0.5],'Rotation',90);
clear xlims ylims;

%set-up subplots for graphs
figureB = figure; set(gcf,'position',[50 400 800 600]);
sub1b = subplot(1,2,1); sub2b = subplot(1,2,2);
figureC = figure; set(gcf,'position',[50 50 800 1000]);

%create dummy matrices to fill with meltwater flux & submerged area for E
%vs W Antarctica in order to fit trendlines to the 2 datasets
WdVdt = []; WAsub = []; EdVdt = []; EAsub = [];
disp('Created figure templates');

%% create subplots & add site info to the map
cd(iceberg_path);
disp('Generating plots...');

%set-up dummy vectors to fill with concatenated variables for all sites
start_yr = []; end_yr = []; avg_x = []; avg_y = []; depth = []; depth_uncert = []; subarea = []; subarea_uncert = [];
meltflux = []; meltflux_uncert = []; meltrate = []; meltrate_uncert = [];
%start & end coordinates for each iceberg
xcoord_o = []; ycoord_o = [];
xcoord_f = []; ycoord_f = [];

%plot
for i = [8:1:length(region) 7:-1:1]
    cd(char(region(i))); disp_names(i) = {strjoin([cellstr(plot_letters(i)),' ',cellstr(leg_names(i))])};
    abbrev_legend(i) = {strjoin([cellstr(leg_abbrevs(i)),' (',cellstr(leg_names(i)),')'])};
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
        decidate_o = convert_to_decimaldate(meltinfo(j).name(end-49:end-42));
        decidate_f = convert_to_decimaldate(meltinfo(j).name(end-34:end-27));
        start_yr = [start_yr; repmat(decidate_o,length(draft),1)]; end_yr = [end_yr; repmat(decidate_f,length(draft),1)];
        plot_yrs = [decidate_o decidate_f];
        clear decidate_*;
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
        plot(Asub,dVdt/86400,plot_marker,'markerfacecolor',region_colors(i,:),'markerfacecolor',region_colors(i,:),...
            'color',region_colors(i,:),'markersize',symbol_size,'linewidth',1); hold on;
        subplot(sub2b);
        errorbar(draft,365*m,365*abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),365*abs(m).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,plot_marker,...
            'markerfacecolor','w','markersize',symbol_size,'markeredgecolor',region_colors(i,:),...
            'color',region_colors(i,:)); hold on;
        %create plot handles for a glacier from each region to plot
        %regional color labels in the legend
        if ~isempty(strmatch('Seller-Bugge',char(region(i))))
            pl(1) = plot(draft,365*m,plot_marker,'markerfacecolor',region_colors(i,:),'markerfacecolor',region_colors(i,:),...
            'color',region_colors(i,:),'markersize',symbol_size,'linewidth',1); hold on;
        elseif ~isempty(strmatch('Thwaites',char(region(i))))
            pl(2) = plot(draft,365*m,plot_marker,'markerfacecolor',region_colors(i,:),'markerfacecolor',region_colors(i,:),...
            'color',region_colors(i,:),'markersize',symbol_size,'linewidth',1); hold on;
        elseif ~isempty(strmatch('Mertz',char(region(i))))
            pl(3) = plot(draft,365*m,plot_marker,'markerfacecolor',region_colors(i,:),'markerfacecolor',region_colors(i,:),...
            'color',region_colors(i,:),'markersize',symbol_size,'linewidth',1); hold on;
        elseif ~isempty(strmatch('Edgeworth-LarsenA',char(region(i))))
            pl(4) = plot(draft,365*m,plot_marker,'markerfacecolor',region_colors(i,:),'markerfacecolor',region_colors(i,:),...
            'color',region_colors(i,:),'markersize',symbol_size,'linewidth',1); hold on;
        else
            plot(draft,365*m,plot_marker,'markerfacecolor',region_colors(i,:),'markerfacecolor',region_colors(i,:),...
            'color',region_colors(i,:),'markersize',symbol_size,'linewidth',1); hold on;
        end
        
        
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
    
    %add symbols to the site maps
    figure(figureA1); colormap(gca,'gray');
    mp1(i) = plot(nanmean(avgx),nanmean(avgy),[map_marker,'k'],'markerfacecolor',region_colors(i,:),'markeredgecolor','k','markersize',12); hold on;
    if strcmp(char(plot_letters(i)),'a)')
        text(nanmean(avgx)-150000,nanmean(avgy)+50000,char(plot_letters(i)),'fontsize',12,'color',[0.5 0.5 0.5]);
    elseif strcmp(char(plot_letters(i)),'b)') || strcmp(char(plot_letters(i)),'c)') || strcmp(char(plot_letters(i)),'e)')
        text(nanmean(avgx)-150000,nanmean(avgy)-50000,char(plot_letters(i)),'fontsize',12,'color',[0.5 0.5 0.5]);
    else
        text(nanmean(avgx)+75000,nanmean(avgy),char(plot_letters(i)),'fontsize',12,'color',[0.5 0.5 0.5]);
    end
    figure(figureA2); colormap(gca,'gray');
    mp2(i) = plot(nanmean(avgx),nanmean(avgy),[map_marker,'k'],'markerfacecolor',region_colors(i,:),'markeredgecolor','k','markersize',12); hold on;
    if strcmp(char(plot_letters(i)),'a)')
        text(nanmean(avgx)-180000,nanmean(avgy)+50000,char(leg_abbrevs(i)),'fontsize',10,'color',[0.5 0.5 0.5]);
    elseif strcmp(char(plot_letters(i)),'b)') || strcmp(char(plot_letters(i)),'c)') || strcmp(char(plot_letters(i)),'e)')
        text(nanmean(avgx)-180000,nanmean(avgy)-50000,char(leg_abbrevs(i)),'fontsize',10,'color',[0.5 0.5 0.5]);
    else
        text(nanmean(avgx)+75000,nanmean(avgy),char(leg_abbrevs(i)),'fontsize',10,'color',[0.5 0.5 0.5]);
    end
       
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
            if strcmp(cellstr(leg_names(i)),'Thwaites')
                set(gca,'xlim',[0 800],'xtick',[0:200:800],'xticklabel',[],...
                    'ylim',[0 70],'ytick',[0:30:70],'yticklabel',[0:30:70],'fontsize',16); grid on;
                figpos = get(gca,'position'); set(gca,'position',[figpos(1) figpos(2)-0.02 figpos(3) figpos(4)]);
            else
                set(gca,'xlim',[0 400],'xtick',[0:100:400],'xticklabel',[],...
                    'ylim',[0 70],'ytick',[0:30:70],'yticklabel',[0:30:70],'fontsize',16); grid on;
            end
        end
    end
    %add axes labels
    if plot_loc(i) == 13
        set(gca,'xlim',[0 400],'xtick',[0:100:400],'xticklabel',[0:100:400]); grid on;
    end
    if plot_loc(i) == 14
        legsub = legend(yrpl,num2str(yrs')); set(legsub,'location','southoutside','NumColumns',5,'fontsize',16);
        legpos = get(legsub,'position'); set(legsub,'position',[0.56 0.10 legpos(3) legpos(4)]); clear legpos;
        set(gca,'xlim',[0 400],'xtick',[0:100:400],'xticklabel',[0:100:400]); grid on;
        xlabel('Draft (m b.s.l.)','fontsize',16);
    end
    if plot_loc(i) == 15
        xlabel('Draft (m b.s.l.)','fontsize',16); ylb = ylabel('Iceberg melt rate (m yr^{-1})','fontsize',16);
        set(gca,'xlim',[0 800],'xtick',[0:200:800],'xticklabel',[0:200:800]);
        set(ylb,'position',[-100 350 -1]);
    end
    text(0.6*max(get(gca,'xlim')),0.85*max(get(gca,'ylim')),[char(plot_letters(i)),' ',char(leg_names(i))],'fontsize',16);
    
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
    saveas(figureD,[figure_path,char(region(i)),'_melt-map.eps'],'epsc'); saveas(figureD,[figure_path,char(region(i)),'_melt-map.png'],'png');
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
%     if i <= 2
%         set(leg,'fontsize',16,'location','southeast');
%     else
        set(leg,'fontsize',16,'location','northwest');
%     end
    saveas(figureE,[figure_path,char(region(i)),'_iceberg-meltrate-v-draft.eps'],'epsc'); saveas(figureE,[figure_path,char(region(i)),'_iceberg-meltrate-v-draft.png'],'png');
    clear xlims ylims;
    
    %save all dates and coords to a tab-delimited text file
    cd ..
    clear subpl indpl daterange;
    close(figureD); close(figureE);
end

%save the overview site maps
figure(figureA1);
[sorted,inds] = sort(leg_ref); mp1_sort = mp1(inds);
legmap = legend(mp1_sort,[char(disp_names(inds))]); set(legmap,'location','northoutside','fontsize',16,'NumColumns',5); 
legmappos = get(legmap,'position'); set(legmap,'position',[0.05 legmappos(2)+0.05 legmappos(3) legmappos(4)]);
gcapos = get(gca,'position'); set(gca,'position',[gcapos(1) 0.09 gcapos(3) gcapos(4)]);
saveas(gcf,[figure_path,'Antarctic-iceberg-map.eps'],'epsc'); saveas(gcf,[figure_path,'Antarctic-iceberg-map.png'],'png');
figure(figureA2);
[sorted,inds] = sort(leg_ref); mp2_sort = mp2(inds);
legmap2 = legend(mp2_sort,[char(abbrev_legend(inds))]); set(legmap2,'location','northoutside','fontsize',12,'NumColumns',5); 
legmap2pos = get(legmap2,'position'); set(legmap2,'position',[0.065 legmap2pos(2)+0.02 legmap2pos(3) legmap2pos(4)]);
gcapos = get(gca,'position'); set(gca,'position',[gcapos(1) 0.09 gcapos(3) gcapos(4)]);
saveas(gcf,[figure_path,'USAPDC-Antarctic-iceberg-map.eps'],'epsc'); saveas(gcf,[figure_path,'USAPDC-Antarctic-iceberg-map.png'],'png');

%save the subplots containing all data
figure(figureB);
subplot(sub1b);
set(gca,'xlim',[0 7.0e6],'xtick',[0:1e6:7e6],'xticklabel',[0:1:7],...
    'ylim',[0 7],'ytick',[0:1:7],'yticklabel',[0:1:7],'fontsize',16); grid on;
xlabel('Submerged area (km^2)','fontsize',16); ylabel('Meltwater flux (m^3 s^{-1})','fontsize',16);
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
text(0.05*max(xlims),0.95*max(ylims),'a) ','color','k','fontsize',16);
subplot(sub2b);
leg1 = legend(pl,[{'West Peninsula (WAP)'};{'West Ice Sheet (WAIS)'};{'East Ice Sheet (EAIS)'};{'East Peninsula (EAP)'}]); 
set(leg1,'location','northeast','fontsize',16,'orientation','vertical'); 
set(gca,'xlim',[0 850],'xtick',[0:200:850],'xticklabel',[0:200:850],...
    'ylim',[0 70],'ytick',[0:10:70],'yticklabel',[0:10:70],'fontsize',16); grid on;
xlabel('Draft (m b.s.l.)','fontsize',16); ylabel('Melt rate (m yr^{-1})','fontsize',16);
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
text(0.05*max(xlims),0.95*max(ylims),'b) ','color','k','fontsize',16);
set(sub1b,'position',[0.1 0.11 0.38 0.750]); set(sub2b,'position',[0.57 0.11 0.38 0.750]);
saveas(gcf,[figure_path,'Antarctic-iceberg-lumped-plots.eps'],'epsc'); saveas(gcf,[figure_path,'Antarctic-iceberg-lumped-plots.png'],'png');

%save the subplots sorted by study site
figure(figureC);
saveas(gcf,[figure_path,'Antarctic-iceberg-subplots.eps'],'epsc'); saveas(gcf,[figure_path,'Antarctic-iceberg-subplots.png'],'png');


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
writetable(T,[iceberg_path,'Antarctic-iceberg-meltinfo.csv']);
disp('Antarctic iceberg melt rate text file written');
clear column_*;
%COORDINATES ONLY
column_names = {'Start Date' 'Start Polar Stereo Easting' 'Start Polar Stereo Northing'...
    'End Date' 'End Polar Stereo Easting' 'End Polar Stereo Northing'};
column_units = {'years' 'meters' 'meters'...
    'years' 'meters' 'meters'};
t=table(start_yr,xcoord_o,ycoord_o,end_yr,xcoord_f,ycoord_f);
t.Properties.VariableNames = column_names; t.Properties.VariableUnits = column_units;
writetable(t,[iceberg_path,'Antarctic-iceberg-PScoords.csv']);
disp('Antarctic iceberg coordinates text file written');
clear T t;
