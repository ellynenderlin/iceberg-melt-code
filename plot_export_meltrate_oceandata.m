%% Plot remotely-sensed iceberg melt vs ocean observations



%% Initialize
clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code','/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean');
% addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/AntarcticMappingTools');

%specify paths & file names for data
CTD_path = '/Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/CTD_Antarctica/';
CTD_data = [CTD_path,'Antarctic-ocean-data.mat'];
RACMO_path = '/Users/ellynenderlin/Research/miscellaneous/RACMO2.3_Antarctica/';
iceberg_path = '/Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt/';
figure_path = [iceberg_path,'figures/'];

%reload compiled data as needed
if ~exist('melt')
    load([iceberg_path,'Antarctic-icebergmelt-comparison.mat']);
end

%specify study site names
region = [{'Edgeworth-LarsenA'},{'Crane-LarsenB'},{'Ronne'},{'Filchner'},{'Amery'},{'Totten'},{'Mertz'},...
    {'Thwaites'},{'Ferrigno-Eltanin'},{'Seller-Bugge'},{'Heim-Marguerite'},{'Widdowson-Biscoe'},{'Cadman-Palmer'},{'Blanchard-Danco'},{'Leonardo-Danco'}];
leg_names = [{'Edgeworth'},{'Crane'},{'Ronne'},{'Filchner'},{'Polar Times'},{'Totten'},{'Mertz'},{'Thwaites'},{'Ferrigno'},{'Seller'},{'Heim'},{'Widdowson'},{'Cadman'},{'Blanchard'},{'Leonardo'}];
for j = 1:length(leg_names); if strcmp(leg_names(j),melt(j).dispname) ~= 1; error('Specified region order must be the same as in the melt structure'); end; end
leg_ref = [9,10,11,12,13,14,15,8,7,6,5,4,3,2,1]; %arrange legend in alphabetical order
region_ref = [1,1,2,2,2,2,2,3,3,4,4,4,4,4,4]; %group according to region (1=EAP, 2=EAIS, 3=WAIS, 4=WAP)

%specify plot params
marker = ['s','s','s','s','s','s','s','s','s','s','s','s','s','s','s']; %modify if you want to change symbols to indicate something about data (E vs W for example)
plot_letters = [{'i)'},{'j)'},{'k)'},{'l)'},{'m)'},{'n)'},{'o)'},{'h)'},{'g)'},{'f)'},{'e)'},{'d)'},{'c)'},{'b)'},{'a)'}]; %plot letters for sites to be used in geographically-arranged subplots
plot_loc = [2,4,6,8,10,12,14,15,13,11,9,7,5,3,1];
EAP_color = [77,172,38]./255; EAIS_color = [184,225,134]./255; WAIS_color = [241,182,218]./255; WAP_color = [208,28,139]./255; 

%create color pallettes
region_colors = zeros(length(region_ref),3);
region_colors(find(region_ref == 1),:) = repmat(EAP_color,length(find(region_ref == 1)),1); 
region_colors(find(region_ref == 2),:) = repmat(EAIS_color,length(find(region_ref == 2)),1); 
region_colors(find(region_ref == 3),:) = repmat(WAIS_color,length(find(region_ref == 3)),1); 
region_colors(find(region_ref == 4),:) = repmat(WAP_color,length(find(region_ref == 4)),1); 
cmap_add = 2.0; cmap_mult = 100; %scalars to convert thermal forcing to symbol color
draft_add = 20; draft_mult = (1/5); %scalars to convert iceberg draft to symbol size
Temp_cmap = cmocean('thermal',2*(cmap_add*cmap_mult)); 
highmelt_cmap = cmocean('amp',100);
Tforcing_cmap = cmocean('amp',600);

%specify generic variables
rho_sw = 1026; %sea water density in kg m^-3
depth_cutoff = 800; %maximum depth of icebergs & therefore ocean observation data of interest in meters
years = [2011.75 2022.25]; year_ticks = [2013:2:2022]; %approximate date range for plots

%specify the buffer region to search for ocean data around the icebergs
buffer = 100000; %m

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

%specify E & W Thwaites and Mertz date indices to fit separate lines to 
%each side of the ice tongues
Thwaites_dates = [melt(strmatch('Thwaites',leg_names)).to melt(strmatch('Thwaites',leg_names)).tf];
[~,~,Thwaites_inds] = unique(Thwaites_dates,'rows');
ThwaitesE = [1,3,5]; ThwaitesW = [2,4,6,7]; %date indices
Thwaites_Erefs = []; Thwaites_Wrefs = [];
for j = ThwaitesE; Thwaites_Erefs = [Thwaites_Erefs; find(Thwaites_inds==j)]; end
for j = ThwaitesW; Thwaites_Wrefs = [Thwaites_Wrefs; find(Thwaites_inds==j)]; end
Mertz_dates = [melt(strmatch('Mertz',leg_names)).to melt(strmatch('Mertz',leg_names)).tf];
[~,~,Mertz_inds] = unique(Mertz_dates,'rows');
MertzE = [4]; MertzW = [1,2,3]; %date indices
Mertz_Erefs = []; Mertz_Wrefs = [];
for j = MertzE; Mertz_Erefs = [Mertz_Erefs; find(Mertz_inds==j)]; end
for j = MertzW; Mertz_Wrefs = [Mertz_Wrefs; find(Mertz_inds==j)]; end
% clear Thwaites_inds ThwaitesE ThwaitesW Mertz_inds MertzE MertzW;

%navigate to the iceberg directory as the default workspace
cd(iceberg_path);
close all;

%% Create profiles of ocean temperature & iceberg depth timeseries (Antarctic-iceberg-oceandata-profiles.eps)

%reload compiled data as needed
if ~exist('melt')
    load([iceberg_path,'Antarctic-icebergmelt-comparison.mat']);
end

%set up plot
close all; 
figure; set(gcf,'position',[50 50 800 1000]);
fprintf('Length of temp colormap = %3i \n',length(Temp_cmap));

%loop through data
disp('Plotting temperature profiles & iceberg depth timeseries');
for i = 1:length(melt)
    if ~isempty(melt(i).oceant)
        subpl = subplot(8,2,plot_loc(i));
        fprintf('Max ocean temp = %4.2f \n',max(max(melt(i).oceanTavg_prof)));
        fprintf('Max temp index = %2i \n',round((max(max(melt(i).oceanTavg_prof))+cmap_add)*cmap_mult));
        
        %plot temperature profiles with colors to distinguish temps
        for j = 1:length(melt(i).oceantavg)
            prof_cmap = []; prof_d = [];
            
            if nanmedian(diff(melt(i).oceandavg_prof)) == 1
                %average over 20 rows to create smoothed profiles at 20 m increments
                for k = 11:20:size(melt(i).oceanTavg_prof(:,j),1)-10
                    prof_d = [prof_d; nanmean(melt(i).oceandavg_prof(k-10:k+10))];
                    if ~isnan(nanmean(melt(i).oceanTavg_prof(k-10:k+10,j))) %&& ~isnan(melt(i).oceanTfreeze_prof(k,j))
                        if round((nanmean(melt(i).oceanTavg_prof(k-10:k+10,j))+cmap_add)*cmap_mult) > 0
                            prof_cmap = [prof_cmap; Temp_cmap(round((nanmean(melt(i).oceanTavg_prof(k-10:k+10,j))+cmap_add)*cmap_mult),:)];
                        else
                            prof_cmap = [prof_cmap; Temp_cmap(1,:)];
                        end
                        %                         plot(melt(i).oceantavg(j),nanmean(melt(i).oceandavg_prof(k-10:k+10)),'.','color',Temp_cmap(round((nanmean(melt(i).oceanTavg_prof(k-10:k+10,j))+3)*100),:),'markersize',10); hold on;
                    else
                        prof_cmap = [prof_cmap; 1,1,1];
                    end
                end
            else
                %plot at native standardized resolution because it should be ~20 m
                for k = 1:size(melt(i).oceanTavg_prof(:,j),1)
                    prof_d = [prof_d; melt(i).oceandavg_prof(k)];
                    if ~isnan(melt(i).oceanTavg_prof(k,j))
                        if round((melt(i).oceanTavg_prof(k,j)+cmap_add)*cmap_mult) > 0
                            prof_cmap = [prof_cmap; Temp_cmap(round((melt(i).oceanTavg_prof(k,j)+cmap_add)*cmap_mult),:)];
                        else
                            prof_cmap = [prof_cmap; Temp_cmap(1,:)];
                        end
                        %                         plot(melt(i).oceantavg(j),melt(i).oceandavg_prof(k),'.','color',Temp_cmap(round((melt(i).oceanTavg_prof(k,j)+3)*100),:),'markersize',10); hold on;
                    else
                        prof_cmap = [prof_cmap; 1,1,1];
                    end
                end
            end
            top_ref = find(sum(prof_cmap,2) < 3,1,'first'); %identify the deepest part of the profile with data as the addition of the plot colors < 3 (sum = 3 is white for NaNs)
            bottom_ref = find(sum(prof_cmap,2) < 3,1,'last'); %identify the deepest part of the profile with data as the addition of the plot colors < 3 (sum = 3 is white for NaNs)
            scatter(repmat(melt(i).oceantavg(j),size(prof_d(top_ref:bottom_ref))),prof_d(top_ref:bottom_ref),14,prof_cmap(top_ref:bottom_ref,:),'filled','s'); hold on;
            clear prof_cmap prof_d bottom_ref;
        end
        
        %plot the iceberg depths for each date
        plot(nanmean([melt(i).to melt(i).tf],2),melt(i).d,'.k'); hold on;
        errorbar(nanmean([melt(i).to melt(i).tf],2),melt(i).d,[],[],abs(nanmean([melt(i).to melt(i).tf],2)-melt(i).to),abs(nanmean([melt(i).to melt(i).tf],2)-melt(i).tf),'.k');
        
        %format plot
        set(gca,'ydir','reverse','xlim',years,'ylim',[0 depth_cutoff],'ytick',[0:250:750],'fontsize',16);
        if plot_loc(i) == 15
            xlabel('Year','fontsize',16); ylabel('Depth (m b.s.l.)','fontsize',16);
        end
        text(min(get(gca,'xlim'))+0.025*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),max(get(gca,'ylim'))-0.2*abs(max(get(gca,'ylim'))-min(get(gca,'ylim'))),[char(plot_letters(i)),' ',char(melt(i).dispname)],'fontsize',16);
        drawnow;

    else
        subpl = subplot(8,2,plot_loc(i));
        %plot the iceberg depths for each date
        plot(nanmean([melt(i).to melt(i).tf],2),melt(i).d,'.k'); hold on;
        errorbar(nanmean([melt(i).to melt(i).tf],2),melt(i).d,[],[],abs(nanmean([melt(i).to melt(i).tf],2)-melt(i).to),abs(nanmean([melt(i).to melt(i).tf],2)-melt(i).tf),'.k');
        
        %format plot
        set(gca,'ydir','reverse','xlim',years,'ylim',[0 depth_cutoff],'ytick',[0:250:750],'fontsize',16);
        text(min(get(gca,'xlim'))+0.025*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),max(get(gca,'ylim'))-0.2*abs(max(get(gca,'ylim'))-min(get(gca,'ylim'))),[char(plot_letters(i)),' ',char(melt(i).dispname)],'fontsize',16);
        drawnow;
        
    end
    
    %format axes
    pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) 1.05*pos(3) 1.15*pos(4)]);
    if plot_loc(i) == 14
        set(gca,'xlim',years,'xtick',year_ticks,'ylim',[0 depth_cutoff],'ytick',[0:250:750],'fontsize',16);
        xlabel('Year','fontsize',16); 
    elseif plot_loc(i) == 15
        set(gca,'xlim',years,'xtick',year_ticks,'ylim',[0 depth_cutoff],'ytick',[0:250:750],'fontsize',16);
        xlabel('Year','fontsize',16); ylbl = ylabel('Depth (m b.s.l.)','fontsize',16);
        set(ylbl,'position',[min(years)-1.5 -3000 -1]);
    else
        set(gca,'xlim',years,'xtick',year_ticks,'xticklabel',[],'ylim',[0 depth_cutoff],'ytick',[0:250:750],'fontsize',16);
    end
    box on;
    
end
%add a colorbar
cbar_start = 0.595; cbar_length = 0.3;
annotation('rectangle',[cbar_start-0.025 0.11 cbar_length+0.05 0.065],'facecolor','w','edgecolor','k');
for j = 1:length(Temp_cmap)
    annotation('line',[cbar_start+j*(cbar_length/length(Temp_cmap)) cbar_start+j*(cbar_length/length(Temp_cmap))],[0.1525 0.170],'color',Temp_cmap(j,:),'linewidth',1.5);
end
annotation('textbox',[cbar_start-0.025 0.135 0.05 0.02],'string',['-',num2str(cmap_add),char(176),'C'],'fontsize',16,'edgecolor','none');
annotation('textbox',[cbar_start-0.025+(0.5*length(Temp_cmap))*(cbar_length/length(Temp_cmap)) 0.135 0.05 0.02],'string',['0',char(176),'C'],'fontsize',16,'edgecolor','none');
annotation('textbox',[cbar_start-0.025+(length(Temp_cmap))*(cbar_length/length(Temp_cmap))-0.005 0.135 0.05 0.02],'string',[num2str(cmap_add),char(176),'C'],'fontsize',16,'edgecolor','none');
annotation('textbox',[cbar_start+0.05 0.115 0.25 0.02],'string','ocean temperature','fontsize',16,'edgecolor','none','fontweight','bold');

%save
saveas(gcf,[figure_path,'Antarctic-iceberg-oceandata-profiles.eps'],'epsc'); saveas(gcf,[figure_path,'Antarctic-iceberg-oceandata-profiles.png'],'png');
disp('iceberg and ocean temp depth profiles saved');

% %create a map that shows the average temperature for each profile down to
% %~100 m depth (Xs) and the median iceberg depth for all sites
% figure(Tm_mapplot);
% for i = 1:length(melt)
%     if ~isempty(melt(i).oceant)
%         for j = 1:length(melt(i).oceantavg)
%             hundred_ref = find(melt(i).oceandavg_prof<=100,1,'last');
%             median_ref = find(melt(i).oceandavg_prof<=nanmedian(melt(i).d),1,'last');
%             if ~isnan(nanmean(melt(i).oceanTavg_prof(1:hundred_ref,j)-melt(i).oceanTfreeze_prof(1:hundred_ref,j)))
%             plot(melt(i).oceanxavg(j),melt(i).oceanyavg(j),'x','color',Temp_cmap(round((nanmean(melt(i).oceanTavg_prof(1:hundred_ref,j)-melt(i).oceanTfreeze_prof(1:hundred_ref,j))+1)*100),:)); hold on;
%             end
%         end
%     end
%     plot(nanmean(melt(i).x),nanmean(melt(i).y),[marker(i),'k'],'markerfacecolor',depth_cmap(round(nanmean(melt(i).d)),:),'markersize',12); hold on;
% end
% %add labels to the location plot
% for i = 1:length(melt)
%     figure(Tm_mapplot); 
%     if strcmp(marker(i),'d')
%         text(nanmean(melt(i).x)+100000,nanmean(melt(i).y),char(plot_letters(i)),'fontsize',16);
%     else
%         if strcmp(char(plot_letters(i)),'f)')
%             text(nanmean(melt(i).x)-200000,nanmean(melt(i).y)-100000,char(plot_letters(i)),'fontsize',16); 
%         elseif strcmp(char(plot_letters(i)),'a)')
%             text(nanmean(melt(i).x)-200000,nanmean(melt(i).y)+100000,char(plot_letters(i)),'fontsize',16); 
%         else
%             text(nanmean(melt(i).x)-200000,nanmean(melt(i).y),char(plot_letters(i)),'fontsize',16); 
%         end
%     end
% end
% %label
% set(gca,'xlim',[-28e5 28e5],'xtick',[-24e5:8e5:24e5],'xticklabel',[-2400:800:2400],...
%     'ylim',[-24e5 24e5],'ytick',[-24e5:8e5:24e5],'yticklabel',[-2400:800:2400],'fontsize',24); grid off;
% xlabel('Easting (km)','fontsize',24); ylabel('Northing (km)','fontsize',24);
% graticuleps(-50:-5:-90,-180:30:180);
% text(0,6.5e5,'85^oS','fontsize',16); text(0,12.0e5,'80^oS','fontsize',16); text(0,17.5e5,'75^oS','fontsize',16); text(0,23.0e5,'70^oS','fontsize',16);
% text(-16.5e5,25.25e5,'-30^oE','fontsize',16); text(12.5e5,25.25e5,'30^oE','fontsize',16); 
% colormap(gca,gray(100001));
% saveas(gcf,'Antarctic-iceberg-oceandata-map.eps','epsc'); saveas(gcf,'Antarctic-iceberg-oceandata-map.png','png');
% %now zoom in on the peninsula and save again
% set(gca,'xlim',[-28e5 -20e5],'xtick',[-28e5:2e5:-20e5],'xticklabel',[-2800:200:-2000],...
%     'ylim',[7.5e5 17.5e5],'ytick',[8e5:2e5:16e5],'yticklabel',[800:200:1600],'fontsize',24);
% graticuleps(-50:-2:-90,-180:10:180);
% saveas(gcf,'AntarcticPeninsula-iceberg-oceandata-map.eps','epsc'); saveas(gcf,'AntarcticPeninsula-iceberg-oceandata-map.png','png');

%% Plot maps and scatterplots of melt rate and ocean thermal forcing
close all; drawnow;
disp('Plotting figures to show iceberg melt rates vs thermal forcing');

%set-up the scatterplots of meltrate vs thermal forcing
Tm_scatterplot = figure; set(gcf,'position',[850 50 800 500]); %observation-based thermal forcing
% TmSOSE_scatterplot = figure; set(gcf,'position',[850 550 800 400]); %observation-based thermal forcing
%knowing the Thwaites data are comprehensive, plot the trendline for those data in the background
figure(Tm_scatterplot); TF_sub = subplot(1,2,1); TFSOSE_sub = subplot(1,2,2);
disp('Thwaites melt vs. thermal forcing:');
%trendlines using observed ocean data
subplot(TF_sub);
%Thwaites East
[TEf,TEgof] = fit(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Erefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Erefs),365*melt(strmatch('Thwaites',leg_names)).m(Thwaites_Erefs),'poly1');
TEci = confint(TEf,0.95); %disp('Thwaites East'); disp(['observed R^2 = ',num2str(round(TEgof.rsquare,2)),' & RMSE = ',num2str(round(TEgof.rmse,1))]);
% plot([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Erefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Erefs))*10)/10],...
%     feval(TEf,[0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Erefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Erefs))*10)/10]),...
%     '-','color',region_colors(strmatch('Thwaites',leg_names),:),'linewidth',2); hold on;
% fill([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Erefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Erefs))*10)/10,...
%     fliplr([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Erefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Erefs))*10)/10])],...
%     [(TEci(1,1).*[0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Erefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Erefs))*10)/10]+TEci(1,2)),...
%     (TEci(2,1).*fliplr([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Erefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Erefs))*10)/10])+TEci(2,2))],...
%     region_colors(strmatch('Thwaites',leg_names),:),'FaceAlpha',0.25,'EdgeColor',region_colors(strmatch('Thwaites',leg_names),:));
%Thwaites West
[TWf,TWgof] = fit(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Wrefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Wrefs),365*melt(strmatch('Thwaites',leg_names)).m(Thwaites_Wrefs),'poly1');
TWci = confint(TWf,0.95); %disp('Thwaites West'); disp(['observed R^2 = ',num2str(round(TWgof.rsquare,2)),' & RMSE = ',num2str(round(TWgof.rmse,1))]);
fprintf('Thwaites W iceberg melt vs observed thermal forcing: %i m/yr per degree (R^2 = %3.2f, RMSE = %3.1f) \n',round(TWf.p1),TWgof.rsquare,TWgof.rmse);
plot([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Wrefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Wrefs))*10)/10],...
    feval(TWf,[0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Wrefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Wrefs))*10)/10]),...
    '-','color',region_colors(strmatch('Thwaites',leg_names),:),'linewidth',2); hold on;
fill([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Wrefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Wrefs))*10)/10,...
    fliplr([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Wrefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Wrefs))*10)/10])],...
    [(TWci(1,1).*[0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Wrefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Wrefs))*10)/10]+TWci(1,2)),...
    (TWci(2,1).*fliplr([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavg(Thwaites_Wrefs)-melt(strmatch('Thwaites',leg_names)).oceanTfp(Thwaites_Wrefs))*10)/10])+TWci(2,2))],...
    region_colors(strmatch('Thwaites',leg_names),:),'FaceAlpha',0.25,'EdgeColor',region_colors(strmatch('Thwaites',leg_names),:));
%trendlines using SOSE data
subplot(TFSOSE_sub);
newer = find(melt(strmatch('Thwaites',leg_names)).to>=2020);
%Thwaites East
[TEf,TEgof] = fit(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Erefs(Thwaites_Erefs < min(newer)))-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Erefs(Thwaites_Erefs < min(newer))),365*melt(strmatch('Thwaites',leg_names)).m(Thwaites_Erefs(Thwaites_Erefs < min(newer))),'poly1');
TEci = confint(TEf,0.95); %disp('Thwaites East'); disp(['modeled R^2 = ',num2str(round(TEgof.rsquare,2)),' & RMSE = ',num2str(round(TEgof.rmse,1))]);
% plot([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Erefs)-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Erefs))*10)/10],...
%     feval(TEf,[0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Erefs)-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Erefs))*10)/10]),...
%     '-','color',region_colors(strmatch('Thwaites',leg_names),:),'linewidth',2); hold on;
% fill([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Erefs)-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Erefs))*10)/10,...
%     fliplr([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Erefs)-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Erefs))*10)/10])],...
%     [(TEci(1,1).*[0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Erefs)-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Erefs))*10)/10]+TEci(1,2)),...
%     (TEci(2,1).*fliplr([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Erefs)-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Erefs))*10)/10])+TEci(2,2))],...
%     region_colors(strmatch('Thwaites',leg_names),:),'FaceAlpha',0.25,'EdgeColor',region_colors(strmatch('Thwaites',leg_names),:));
%Thwaites West
[TWf,TWgof] = fit(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Wrefs(Thwaites_Wrefs < min(newer)))-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Wrefs(Thwaites_Wrefs < min(newer))),365*melt(strmatch('Thwaites',leg_names)).m(Thwaites_Wrefs(Thwaites_Wrefs < min(newer))),'poly1');
TWci = confint(TWf,0.95); %disp('Thwaites West'); disp(['modeled R^2 = ',num2str(round(TWgof.rsquare,2)),' & RMSE = ',num2str(round(TWgof.rmse,1))]);
% plot([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Wrefs)-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Wrefs))*10)/10],...
%     feval(TWf,[0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Wrefs)-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Wrefs))*10)/10]),...
%     '-','color',region_colors(strmatch('Thwaites',leg_names),:),'linewidth',2); hold on;
% fill([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Wrefs)-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Wrefs))*10)/10,...
%     fliplr([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Wrefs)-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Wrefs))*10)/10])],...
%     [(TWci(1,1).*[0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Wrefs)-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Wrefs))*10)/10]+TWci(1,2)),...
%     (TWci(2,1).*fliplr([0:0.1:ceil(max(melt(strmatch('Thwaites',leg_names)).oceanTavgSOSE(Thwaites_Wrefs)-melt(strmatch('Thwaites',leg_names)).oceanTfpSOSE(Thwaites_Wrefs))*10)/10])+TWci(2,2))],...
%     region_colors(strmatch('Thwaites',leg_names),:),'FaceAlpha',0.25,'EdgeColor',region_colors(strmatch('Thwaites',leg_names),:));
clear newer;

%set-up full Antarctic map
Tm_mapplot = figure; set(gcf,'position',[50 50 800 800]);
im_cmap = colormap(gray(10001)); im_cmap(1,:) = [1 1 1];
imagesc(IM.x,IM.y,IM.z); hold on; axis xy equal; colormap(gca,im_cmap);
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

%set-up the map figures for 4 regions with ocean data: Antarctic Peninsula, Thwaites, Mertz, Filchner
%AP map
Tm_mapplot_AP = figure; set(gcf,'position',[50 50 800 800]);
imagesc(IM.x,IM.y,IM.z); hold on; axis xy equal; colormap(gca,im_cmap);
xlims = [-2.58e6 -2.08e6]; ylims = [8.5e5 16.0e5];
set(gca,'xlim',xlims,'ylim',ylims,...
    'xtick',[(ceil(min(xlims)/1000)*1000):50000:(floor(max(xlims)/1000)*1000)],...
    'xticklabel',[(ceil(min(xlims)/1000)):50:(floor(max(xlims)/1000))],...
    'ytick',[(ceil(min(ylims)/1000)*1000):50000:(floor(max(ylims)/1000)*1000)],...
    'yticklabel',[(ceil(min(ylims)/1000)):50:(floor(max(ylims)/1000))],...
    'XTickLabelRotation',45,'fontsize',16); %grid on;
xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
figpos = get(gcf,'position'); set(gcf,'position',[figpos(1) figpos(2) figpos(3) range(ylims)/range(xlims)*figpos(3)+25]); %resize vertical dimension to maximize figure window usage
graticuleps(-60:-2:-70,-90:5:-30);
text(-24.15e5,max(ylims)-3e4,['64',char(176),'S'],'fontsize',16,'Rotation',55,'Color',[0.5 0.5 0.5]); 
text(-21.375e5,max(ylims)-3e4,['66',char(176),'S'],'fontsize',16,'Rotation',55,'Color',[0.5 0.5 0.5]); 
text(max(xlims)-2.5e4,12.1e5,['68',char(176),'S'],'fontsize',16,'Rotation',55,'Color',[0.5 0.5 0.5]);
text(-22.68e5,max(ylims)-1e3,['-55',char(176),'E'],'fontsize',16,'Rotation',-35,'Color',[0.5 0.5 0.5]); 
text(min(xlims)+3e3,14.99e5,['-60',char(176),'E'],'fontsize',16,'Rotation',-30,'Color',[0.5 0.5 0.5]); 
text(min(xlims)+3e3,12.115e5,['-65',char(176),'E'],'fontsize',16,'Rotation',-25.0,'Color',[0.5 0.5 0.5]);
text(min(xlims)+3e3,9.46e5,['-70',char(176),'E'],'fontsize',16,'Rotation',-20,'Color',[0.5 0.5 0.5]);
clear xlims ylims; 
%Thwaites Glacier map
Tm_mapplot_TG = figure; set(gcf,'position',[250 50 800 800]);
imagesc(IM.x,IM.y,IM.z); hold on; axis xy equal; colormap(gca,im_cmap);
xlims = [-1.7e6 -1.525e6]; ylims = [-6e5 -3.25e5];
set(gca,'xlim',xlims,'ylim',ylims,...
    'xtick',[(ceil(min(xlims)/1000)*1000):50000:(floor(max(xlims)/1000)*1000)],...
    'xticklabel',[(ceil(min(xlims)/1000)):50:(floor(max(xlims)/1000))],...
    'ytick',[(ceil(min(ylims)/1000)*1000):50000:(floor(max(ylims)/1000)*1000)],...
    'yticklabel',[(ceil(min(ylims)/1000)):50:(floor(max(ylims)/1000))],...
    'XTickLabelRotation',45,'fontsize',16); %grid on;
xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
figpos = get(gcf,'position'); set(gcf,'position',[figpos(1) figpos(2) figpos(3) range(ylims)/range(xlims)*figpos(3)+25]); %resize vertical dimension to maximize figure window usage
graticuleps(-72:-1:-78,-110:2:-100);
text(-16.5e5,min(ylims)+1.25e4,['74',char(176),'S'],'fontsize',16,'Rotation',-72.5,'Color',[0.5 0.5 0.5]); 
text(-15.32e5,min(ylims)+1.25e4,['75',char(176),'S'],'fontsize',16,'Rotation',-70,'Color',[0.5 0.5 0.5]); 
text(max(xlims)-1.5e4,max(ylims)-0.45e4,['-102',char(176),'E'],'fontsize',16,'Rotation',10,'Color',[0.5 0.5 0.5]); 
text(max(xlims)-1.5e4,-3.865e5,['-104',char(176),'E'],'fontsize',16,'Rotation',12.5,'Color',[0.5 0.5 0.5]); 
text(max(xlims)-1.5e4,-4.44e5,['-106',char(176),'E'],'fontsize',16,'Rotation',15,'Color',[0.5 0.5 0.5]);
text(max(xlims)-1.5e4,-5.025e5,['-108',char(176),'E'],'fontsize',16,'Rotation',17.5,'Color',[0.5 0.5 0.5]);
text(max(xlims)-1.5e4,-5.630e5,['-110',char(176),'E'],'fontsize',16,'Rotation',20,'Color',[0.5 0.5 0.5]);
clear xlims ylims; 
%Mertz Ice Tongue map
Tm_mapplot_MI = figure; set(gcf,'position',[450 50 800 800]);
imagesc(IM.x,IM.y,IM.z); hold on; axis xy equal; colormap(gca,im_cmap);
xlims = [1.25e6 1.5e6]; ylims = [-21.5e5 -20.0e5];
set(gca,'xlim',xlims,'ylim',ylims,...
    'xtick',[(ceil(min(xlims)/1000)*1000):50000:(floor(max(xlims)/1000)*1000)],...
    'xticklabel',[(ceil(min(xlims)/1000)):50:(floor(max(xlims)/1000))],...
    'ytick',[(ceil(min(ylims)/1000)*1000):50000:(floor(max(ylims)/1000)*1000)],...
    'yticklabel',[(ceil(min(ylims)/1000)):50:(floor(max(ylims)/1000))],...
    'XTickLabelRotation',45,'fontsize',16); %grid on;
xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
figpos = get(gcf,'position'); set(gcf,'position',[figpos(1) figpos(2) figpos(3) range(ylims)/range(xlims)*figpos(3)+25]); %resize vertical dimension to maximize figure window usage
graticuleps(-66:-1:-70,142:2:150);
text(min(xlims)+0.75e4,max(ylims)-0.45e4,['148',char(176),'E'],'fontsize',16,'Rotation',-58,'Color',[0.5 0.5 0.5]); 
text(13.57e5,max(ylims)-0.45e4,['146',char(176),'E'],'fontsize',16,'Rotation',-56,'Color',[0.5 0.5 0.5]);
text(14.625e5,max(ylims)-0.45e4,['144',char(176),'E'],'fontsize',16,'Rotation',-54,'Color',[0.5 0.5 0.5]);
text(13.3e5,max(ylims)-1.5e4,['68',char(176),'S'],'fontsize',16,'Rotation',30,'Color',[0.5 0.5 0.5]);
text(14.85e5,-20.45e5,['67',char(176),'S'],'fontsize',16,'Rotation',32,'Color',[0.5 0.5 0.5]);
clear xlims ylims; 
%Totten Ice Shelf map
Tm_mapplot_TI = figure; set(gcf,'position',[450 50 800 800]);
imagesc(IM.x,IM.y,IM.z); hold on; axis xy equal; colormap(gca,im_cmap);
xlims = [2.20e6 2.45e6]; ylims = [-12.5e5 -11.0e5];
set(gca,'xlim',xlims,'ylim',ylims,...
    'xtick',[(ceil(min(xlims)/1000)*1000):50000:(floor(max(xlims)/1000)*1000)],...
    'xticklabel',[(ceil(min(xlims)/1000)):50:(floor(max(xlims)/1000))],...
    'ytick',[(ceil(min(ylims)/1000)*1000):50000:(floor(max(ylims)/1000)*1000)],...
    'yticklabel',[(ceil(min(ylims)/1000)):50:(floor(max(ylims)/1000))],...
    'XTickLabelRotation',45,'fontsize',16); %grid on;
xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
figpos = get(gcf,'position'); set(gcf,'position',[figpos(1) figpos(2) figpos(3) range(ylims)/range(xlims)*figpos(3)+25]); %resize vertical dimension to maximize figure window usage
graticuleps(-64:-1:-68,112:2:120);
text(2.27e6,max(ylims)-0.35e4,['116',char(176),'E'],'fontsize',16,'Rotation',-24,'Color',[0.5 0.5 0.5]); 
text(min(xlims)+0.3e4,-1.1675e6,['118',char(176),'E'],'fontsize',16,'Rotation',-28,'Color',[0.5 0.5 0.5]);
text(2.208e6,min(ylims)+3e3,['67',char(176),'S'],'fontsize',16,'Rotation',60,'Color',[0.5 0.5 0.5]);
text(2.3375e6,min(ylims)+3e3,['66',char(176),'S'],'fontsize',16,'Rotation',60,'Color',[0.5 0.5 0.5]);
clear xlims ylims; 
%Filchner Ice Shelf map
Tm_mapplot_FI = figure; set(gcf,'position',[650 50 800 800]);
imagesc(IM.x,IM.y,IM.z); hold on; axis xy equal; colormap(gca,im_cmap);
xlims = [-0.875e6 -0.725e6]; ylims = [9.75e5 11.75e5];
set(gca,'xlim',xlims,'ylim',ylims,...
    'xtick',[(ceil(min(xlims)/1000)*1000):50000:(floor(max(xlims)/1000)*1000)],...
    'xticklabel',[(ceil(min(xlims)/1000)):50:(floor(max(xlims)/1000))],...
    'ytick',[(ceil(min(ylims)/1000)*1000):50000:(floor(max(ylims)/1000)*1000)],...
    'yticklabel',[(ceil(min(ylims)/1000)):50:(floor(max(ylims)/1000))],...
    'XTickLabelRotation',45,'fontsize',16); %grid on;
xlabel('Easting (km)','fontsize',16); ylabel('Northing (km)','fontsize',16);
figpos = get(gcf,'position');set(gcf,'position',[figpos(1) figpos(2) figpos(3) range(ylims)/range(xlims)*figpos(3)+25]); %resize vertical dimension to maximize figure window usage
graticuleps(-76:-1:-80,-42:2:-32);
text(max(xlims)-0.85e4,10.825e5,['-34',char(176),'E'],'fontsize',16,'Rotation',-60,'Color',[0.5 0.5 0.5]); 
text(-7.705e5,10.570e5,['-36',char(176),'E'],'fontsize',16,'Rotation',-56,'Color',[0.5 0.5 0.5]);
text(-8.070e5,10.290e5,['-38',char(176),'E'],'fontsize',16,'Rotation',-52,'Color',[0.5 0.5 0.5]);
text(-8.425e5,10.0e5,['-40',char(176),'E'],'fontsize',16,'Rotation',-48,'Color',[0.5 0.5 0.5]);
text(-8.075e5,max(ylims)-0.65e4,['78',char(176),'S'],'fontsize',16,'Rotation',32,'Color',[0.5 0.5 0.5]);
text(max(xlims)-1.0e4,10.85e5,['77',char(176),'S'],'fontsize',16,'Rotation',34,'Color',[0.5 0.5 0.5]);
clear xlims ylims; 

%set up a colormap for iceberg draft
depth_cmap_top = cmocean('-topo',round(depth_cutoff/2)); 
depth_cmap_bottom = cmocean('-topo',round(depth_cutoff*1.5)); 
depth_cmap = [depth_cmap_top(1:floor(depth_cutoff/4),:); depth_cmap_bottom(ceil(depth_cutoff*0.75):end,:)];
clear depth_cmap_*;
%depth_cmap = cmocean('deep',depth_cutoff); 

tempref = [];
for i = 1:length(melt)
    disp(char(leg_names(i)));
    disp_names(i) = {strjoin([cellstr(plot_letters(i)),' ',cellstr(leg_names(i))])};
    landsats = dir([iceberg_path,char(leg_names(i)),'/LC*PS.TIF']);
    sitex = []; sitey = []; %set up empty cells to insert ocean data coordinates (if data exist) & iceberg coordinates for adjusting site maps

    %identify the maximum iceberg draft
    for k = 1:length(melt(i).x)
        draft_map(k,:) = depth_cmap(round(nanmean(melt(i).d(k))),:);
        symbol_color(k) = ceil(2*(nanmedian(melt(i).m(k))*365)); 
        if symbol_color(k) > length(highmelt_cmap); symbol_color(k) = length(highmelt_cmap); end
    end
    max_ind = find(melt(i).d == max(melt(i).d)); %identify the deepest iceberg
    
    %find appropriate ocean data & extract thermal forcing estimates
    if ~isempty(melt(i).oceant)
        TFcoords = []; TF = [];
        
        %loop through remotely-sensed data and extract ocean forcing information for each iceberg
        for j = 1:length(melt(i).to) %length of melt(i).to corresponds to the number of icebergs
            %identify the time span of remotely-sensed iceberg melt rate estimates
            %(bi-annual=2, annual=1, or seasonal=0)
            if melt(i).tf(j)-melt(i).to(j) >= 2
                timespan = 2;
            elseif melt(i).tf(j)-melt(i).to(j) >= 1
                timespan = 1;
            else
                timespan = 0;
            end
            
            %if seasonal, find ocean data from approximately the same season
            %minrefs specifies the indices for the closest date (there may be multiple profiles) & oceantemps and oceansals are the corresponding profiles
            if size(melt(i).oceanx,2) ~=1; melt(i).oceanx = melt(i).oceanx'; melt(i).oceany = melt(i).oceany'; end %make sure coordinates are always a column vector
            if timespan == 0
                deciseas = nanmean([melt(i).to(j) melt(i).tf(j)]-floor(melt(i).to(j)),2); if deciseas > 1; deciseas = deciseas - floor(deciseas); end
                
                [mindiff,minref] = min(abs((melt(i).oceant-floor(melt(i).oceant))-deciseas));
                if melt(i).oceant(minref)-floor(melt(i).oceant(minref)) >= melt(i).to(j)-floor(melt(i).to(j)) && melt(i).oceant(minref)-floor(melt(i).oceant(minref)) <= melt(i).tf(j)-floor(melt(i).to(j)) %if to and tf are in the same year & minref is in between, find all between
                    minrefs = find(melt(i).oceant-floor(melt(i).oceant(minref)) >= melt(i).to(j)-floor(melt(i).to(j)) & melt(i).oceant-floor(melt(i).oceant(minref)) <= melt(i).tf(j)-floor(melt(i).to(j)));
                elseif melt(i).oceant(minref)-floor(melt(i).oceant(minref)) <= melt(i).to(j)-floor(melt(i).to(j)) && melt(i).oceant(minref)-floor(melt(i).oceant(minref)) <= melt(i).tf(j)-floor(melt(i).tf(j)) %if tf is in a different year than to & minref is in between, find all between
                    minrefs = find(melt(i).oceant-floor(melt(i).oceant(minref)) <= melt(i).to(j)-floor(melt(i).to(j)) & melt(i).oceant-floor(melt(i).oceant(minref)) <= melt(i).tf(j)-floor(melt(i).tf(j)));
                else
                    if mindiff < 0.5 %if there are no data that fall within the seasonal range of to and tf, find data within +/-3 months of the central day of year
                        minrefs = find(abs((melt(i).oceant-floor(melt(i).oceant))-deciseas) <= 4/12);
                    else
                        minrefs = find(abs((melt(i).oceant-floor(melt(i).oceant))-deciseas) <= mindiff + 1/12);
                    end
                end
                oceanx = melt(i).oceanx(minrefs); oceany = melt(i).oceany(minrefs);
                oceantemps = melt(i).oceanT(:,minrefs); oceansals = melt(i).oceanS(:,minrefs); oceandepths = melt(i).oceand(:,minrefs);
                clear minref deciseas mindiff;
            else
                %if annual or bi-annual, find the closest year of ocean data
                [~,minref] = min(abs(melt(i).oceant-nanmean([melt(i).to(j) melt(i).tf(j)])));
                if melt(i).oceant(minref)-nanmean([melt(i).to(j) melt(i).tf(j)]) > 0
                    minrefs = find(melt(i).oceant>=melt(i).oceant(minref) & melt(i).oceant<=melt(i).oceant(minref)+timespan);
                else
                    minrefs = find(melt(i).oceant<=melt(i).oceant(minref) & melt(i).oceant>=melt(i).oceant(minref)-timespan);
                end

                oceanx = melt(i).oceanx(minrefs); oceany = melt(i).oceany(minrefs);
                oceantemps = melt(i).oceanT(:,minrefs); oceansals = melt(i).oceanS(:,minrefs); oceandepths = melt(i).oceand(:,minrefs);
                clear minref;
                
            end
            
            %extract temperature metrics over the iceberg draft
            for k = 1:length(minrefs)
                %identify the bottom of each profile
                if ~isempty(find(oceandepths(:,k)<=melt(i).d(j),1,'last'))
                    bottomrefs(k) = find(oceandepths(:,k)<=melt(i).d(j),1,'last'); %index for the deepest observation for each profile
                    bottomT(k) = melt(i).oceanT(bottomrefs(k),minrefs(k)); bottomS(k) = melt(i).oceanS(bottomrefs(k),minrefs(k));
                else
                    bottomrefs(k) = NaN; bottomT(k) = NaN; bottomS(k) = NaN;
                end
                
                %use the trapz function to calculate the mean for each profile
                %if its maximum observation depth is >90% of the iceberg draft
                if 0.9*max(oceandepths(~isnan(oceantemps(:,k)),k)) > melt(i).d(j) & min(oceandepths(~isnan(oceantemps(:,k)),k)) < 50
                    Tavg(k) = vertmean2(-oceandepths(:,k),oceantemps(:,k),-melt(i).d(j));
                else
                    Tavg(k) = NaN;
                end
                %repeat averaging but for salinity (may not have salinity corresponding to all temp observations)
                if 0.9*max(oceandepths(~isnan(oceansals(:,k)),k)) > melt(i).d(j) & min(oceandepths(~isnan(oceansals(:,k)),k)) < 50
                    Savg(k) = vertmean2(-oceandepths(:,k),oceansals(:,k),-melt(i).d(j));
                else
                    Savg(k) = NaN;
                end
                %repeat averaging, but to calculate the average freezing temperature of sea water (Tfp = -5.73*10^-2 (C/psu)*salinity + 8.32*10^-2 (C) - 7.61*10^-4 (C/dbar)*pressure)
                %pressure is approximately equivalent to depth
                if 0.9*max(oceandepths(~isnan(oceansals(:,k)),k)) > melt(i).d(j) & min(oceandepths(~isnan(oceansals(:,k)),k)) < 50
                    Tfpavg(k) = vertmean2(-oceandepths(:,k),(((-5.73*10^-2).*oceansals(:,k)) + (8.32*10^-2) - ((7.61*10^-4).*oceandepths(:,k))),-melt(i).d(j));
                else
                    Tfpavg(k) = NaN;
                end
            end

            
            %set the size of the melt rate vs thermal forcing scatterplot
            %symbols so that they vary with draft
            draft_size(j) = round(melt(i).d(j).*draft_mult)+draft_add;
            
            %compile thermal forcing data to plot only one thermal forcing 
            %estimate for all icebergs that use the same profile
            TFcoords = [TFcoords; oceanx oceany]; TF = [TF; (Tavg - Tfpavg)'];
            
            clear minrefs oceantemps oceansals oceandepths oceanx oceany timespan Tavg Tfpavg Savg;
        end
        
        %plot the median thermal forcing from each profile on the overview map
        figure(Tm_mapplot);
        %MEDIAN OF AVERAGE THERMAL FORCING FOR ALL ICEBERGS
        %         %calculate the median thermal forcing for each profile
        %         [unique_coords,unique_refs,inds] = unique(TFcoords,'rows');
        %         for j = 1:max(unique_refs)
        %             %median of thermal forcing
        %             TFmedian(j) = nanmedian(TF(inds==j));
        %             
        %             %create the colormap to show thermal forcing on the map
        %             if ~isnan(TFmedian(j))
        %                 tempref = round(TFmedian(j)*(2*cmap_mult));
        %                 if tempref < 1
        %                     temp_map(j,:) = [0 0 0];
        %                 else
        %                     temp_map(j,:) = Tforcing_cmap(tempref,:);
        %                 end
        %                 clear tempref;
        %             else
        %                 temp_map(j,:) = [1 1 1];
        %             end
        %         end
        %         sitex = [sitex; unique_coords(sum(temp_map,2)~=3,1)]; sitey = [sitey; unique_coords(sum(temp_map,2)~=3,2)]; 
        %         scatter(unique_coords(sum(temp_map,2)~=3,1),unique_coords(sum(temp_map,2)~=3,2),16,temp_map((sum(temp_map,2)~=3),:),'filled','o'); hold on;
        %MEDIAN THERMAL FORCING OVER MAXIMUM ICEBERG DRAFT FOR MAP VISUALIZATION ONLY
        for k = 1:size(melt(i).oceanT,2)
            %calculate the average temperature for the profile
            if 0.9*max(melt(i).oceand(~isnan(melt(i).oceanT(:,k)),k)) > melt(i).d(max_ind) & min(melt(i).oceand(~isnan(melt(i).oceanT(:,k)),k)) < 50
                Tavg(k) = vertmean2(-melt(i).oceand(:,k),melt(i).oceanT(:,k),-melt(i).d(max_ind));
            else
                Tavg(k) = NaN;
            end
            %calculate the average freezing temperature for the profile
            if 0.9*max(melt(i).oceand(~isnan(melt(i).oceanS(:,k)),k)) > melt(i).d(max_ind) & min(melt(i).oceand(~isnan(melt(i).oceanS(:,k)),k)) < 50
                Tfpavg(k) = vertmean2(-melt(i).oceand(:,k),(((-5.73*10^-2).*melt(i).oceanS(:,k)) + (8.32*10^-2) - ((7.61*10^-4).*melt(i).oceand(:,k))),-melt(i).d(max_ind));
            else
                Tfpavg(k) = NaN;
            end
            Tfavg(k) = Tavg(k) - Tfpavg(k); %thermal forcing
            
            %assign to a color for scatterplot
            if ~isnan(Tfavg(k))
%                 tempref = round(Tfavg(k)*(2*cmap_mult)); %temp index if using Tforcing_cmap
                tempref = round((Tavg(k)+cmap_add)*cmap_mult); %temp index if using Temp_cmap
                if tempref < 1
                    temp_map(k,:) = [0 0 0];
                else
%                     temp_map(k,:) = Tforcing_cmap(tempref,:);
                    temp_map(k,:) = Temp_cmap(tempref,:);
                end
                clear tempref;
            else
                temp_map(k,:) = [1 1 1];
            end
        end
        scatter(melt(i).oceanx(sum(temp_map,2)~=3),melt(i).oceany(sum(temp_map,2)~=3),16,temp_map((sum(temp_map,2)~=3),:),'filled','o'); hold on;
        clear Tavg Tfpavg Tfavg;
        
        %plot data on the regional map
        if region_ref(i) == 1 || region_ref(i) == 4
            figure(Tm_mapplot_AP); %call the AP map figure handle
        elseif strcmp(cellstr(leg_names(i)),'Thwaites')
            figure(Tm_mapplot_TG); %call the Thwaites map figure handle
        elseif ~isempty(strmatch('Mertz',char(leg_names(i))))
            figure(Tm_mapplot_MI); %call the Mertz map figure handle
        elseif ~isempty(strmatch('Totten',char(leg_names(i))))
            figure(Tm_mapplot_TI); %call the Totten map figure handle
        elseif ~isempty(strmatch('Filchner',char(leg_names(i))))
            figure(Tm_mapplot_FI); %call the Filchner map figure handle
        end
        scatter(melt(i).oceanx(sum(temp_map,2)~=3),melt(i).oceany(sum(temp_map,2)~=3),20,temp_map((sum(temp_map,2)~=3),:),'filled','o','markeredgecolor','w'); hold on; %ocean temperature over maximum draft
%         scatter(melt(i).x,melt(i).y,round(melt(i).d/5)+20,highmelt_cmap(symbol_color,:),'filled','s','markeredgecolor','k','linewidth',0.5); hold on; %icebergs color-coded by meltrate
        scatter(melt(i).x,melt(i).y,round(melt(i).d/5)+20,region_colors(region_ref(i),:),'filled','s','markeredgecolor','k','linewidth',0.5,'markerfacealpha',1); hold on; %icebergs color-coded by region
        drawnow; clear temp_map;
        
        
        %create an individual plot of thermal forcing vs meltrate for the study site
        site_scatter = figure; set(gcf,'position',[850 650 800 400]);
        scatter(melt(i).oceanTavg-melt(i).oceanTfp,100*melt(i).m,2*draft_size,region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
        set(gca,'fontsize',16); grid on;
        xlabel(['Thermal forcing (',char(176),'C above freezing)'],'fontsize',16); ylabel('Melt rate (cm/d)','fontsize',16);
        title(melt(i).dispname); drawnow;
        [~,gof] = fit(melt(i).oceanTavg-melt(i).oceanTfp,melt(i).m,'poly1'); %disp(['Trendline r^2 = ',num2str(gof.rsquare)]);
        if gof.rsquare > 0.5
            saveas(gcf,[figure_path,char(melt(i).name),'-iceberg-oceandata-scatterplot.eps'],'epsc'); saveas(gcf,[figure_path,char(melt(i).name),'-iceberg-oceandata-scatterplot.png'],'png');
        end
        close(site_scatter); clear gof;
        
        %add observation-based thermal forcing vs meltrate data to composite scatterplot
        figure(Tm_scatterplot); subplot(TF_sub)
        if ~isempty(strmatch('Edgeworth',char(leg_names(i)))) %only create a handle for one study site from each region (EAP,EAIS,WAIS,WAP)
            sp(1) = scatter(melt(i).oceanTavg-melt(i).oceanTfp,365*melt(i).m,draft_size,region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
        elseif ~isempty(strmatch('Mertz',char(leg_names(i))))
            sp(2) = scatter(melt(i).oceanTavg-melt(i).oceanTfp,365*melt(i).m,draft_size,region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
        elseif ~isempty(strmatch('Thwaites',char(leg_names(i))))
            sp(3) = scatter(melt(i).oceanTavg-melt(i).oceanTfp,365*melt(i).m,draft_size,region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
        elseif ~isempty(strmatch('Cadman',char(leg_names(i))))
            sp(4) = scatter(melt(i).oceanTavg-melt(i).oceanTfp,365*melt(i).m,draft_size,region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
        else
            scatter(melt(i).oceanTavg-melt(i).oceanTfp,365*melt(i).m,draft_size,region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
        end
        %         for j = 1:length(melt(i).to)
        %             plot(melt(i).oceanTavg(j)-melt(i).oceanTfp(j,1),365*melt(i).m(j),[marker(i),'k'],'markerfacecolor',depth_cmap(round(melt(i).d(j)),:),'markersize',12); hold on;
        %         end
        
    else
        %set symbol sizes so that they vary with draft
        for j = 1:length(melt(i).to)
            draft_size(j) = round(melt(i).d(j)*draft_mult)+draft_add;
        end
        
        %plot data on the regional map
        if region_ref(i) == 1 || region_ref(i) == 4
            figure(Tm_mapplot_AP); %call the AP map figure handle
        elseif strcmp(cellstr(leg_names(i)),'Thwaites')
            figure(Tm_mapplot_TG); %call the Thwaites map figure handle
        elseif ~isempty(strmatch('Mertz',char(leg_names(i))))
            figure(Tm_mapplot_MI); %call the Mertz map figure handle
        elseif ~isempty(strmatch('Totten',char(leg_names(i))))
            figure(Tm_mapplot_TI); %call the Totten map figure handle
        elseif ~isempty(strmatch('Filchner',char(leg_names(i))))
            figure(Tm_mapplot_FI); %call the Filchner map figure handle
        end
%         scatter(melt(i).x,melt(i).y,round(melt(i).d/5)+20,highmelt_cmap(symbol_color,:),'filled','s','markeredgecolor','k','linewidth',0.5); hold on; %icebergs color-coded by meltrate
        scatter(melt(i).x,melt(i).y,round(melt(i).d/5)+20,region_colors(region_ref(i),:),'filled','s','markeredgecolor','k','linewidth',0.5,'markerfacealpha',1); hold on; %icebergs color-coded by region
        drawnow;
    end

    %add SOSE modeled-thermal forcing vs meltrate to scatterplot
%     figure(TmSOSE_scatterplot);
    older = find(melt(i).to<2023);
    figure(Tm_scatterplot); subplot(TFSOSE_sub);
    if ~isempty(strmatch('Edgeworth',char(leg_names(i)))) %only create a handle for one study site from each region (EAP,EAIS,WAIS,WAP)
        spS(1) = scatter(melt(i).oceanTavgSOSE(older)-melt(i).oceanTfpSOSE(older),365*melt(i).m(older),draft_size(older),region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
    elseif ~isempty(strmatch('Mertz',char(leg_names(i))))
        spS(2) = scatter(melt(i).oceanTavgSOSE(older)-melt(i).oceanTfpSOSE(older),365*melt(i).m(older),draft_size(older),region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
    elseif ~isempty(strmatch('Thwaites',char(leg_names(i))))
        spS(3) = scatter(melt(i).oceanTavgSOSE(older)-melt(i).oceanTfpSOSE(older),365*melt(i).m(older),draft_size(older),region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
    elseif ~isempty(strmatch('Cadman',char(leg_names(i))))
        spS(4) = scatter(melt(i).oceanTavgSOSE(older)-melt(i).oceanTfpSOSE(older),365*melt(i).m(older),draft_size(older),region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
    else
        scatter(melt(i).oceanTavgSOSE(older)-melt(i).oceanTfpSOSE(older),365*melt(i).m(older),draft_size(older),region_colors(i,:),'filled','s','markeredgecolor','k'); hold on;
    end
    clear draft_size; clear older;
    
    %add data to the overview map
    figure(Tm_mapplot);
    mp(i) = scatter(melt(i).x(max_ind),melt(i).y(max_ind),round(melt(i).d(max_ind)/5)+10,highmelt_cmap(symbol_color(max_ind),:),'filled','s','markeredgecolor','k','linewidth',0.5); hold on; %only plot a symbol for the deepest-drafted iceberg
    clear temp_map draft_map symbol_color;

end

%format the overview map & save
figure(Tm_mapplot);
%add labels
for i = 1:length(melt)
    if strcmp(char(plot_letters(i)),'a)')
        text(nanmean(melt(i).x)-150000,nanmean(melt(i).y)+50000,char(plot_letters(i)),'fontsize',12,'color',[0.5 0.5 0.5]);
    elseif strcmp(char(plot_letters(i)),'b)') || strcmp(char(plot_letters(i)),'c)') || strcmp(char(plot_letters(i)),'e)')
        text(nanmean(melt(i).x)-150000,nanmean(melt(i).y)-50000,char(plot_letters(i)),'fontsize',12,'color',[0.5 0.5 0.5]);
    else
        text(nanmean(melt(i).x)+75000,nanmean(melt(i).y),char(plot_letters(i)),'fontsize',12,'color',[0.5 0.5 0.5]);
    end
end
%add color & size legends for iceberg data
rectangle('position',[-26.0e5 -22.0e5 28e5 15e5],'facecolor','w','edgecolor','k'); xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
%iceberg sizes
scatter(min(xlims)+0.4*range(xlims),min(ylims)+0.07*range(ylims),round(100*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(min(xlims)+0.425*range(xlims),min(ylims)+0.07*range(ylims),'100 m','fontsize',16); %scaling y-offset = 0.16*range(ylims)
scatter(min(xlims)+0.4*range(xlims),min(ylims)+0.15*range(ylims),round(300*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(min(xlims)+0.425*range(xlims),min(ylims)+0.15*range(ylims),'300 m','fontsize',16); %scaling y-offset = 0.11*range(ylims)
scatter(min(xlims)+0.4*range(xlims),min(ylims)+0.23*range(ylims),round(500*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(min(xlims)+0.425*range(xlims),min(ylims)+0.98*0.23*range(ylims),'500 m','fontsize',16); %scaling y-offset = 0.05*range(ylims)
text(min(xlims)+0.40*range(xlims),min(ylims)+0.305*range(ylims),'iceberg','fontsize',16,'fontweight','bold');
text(min(xlims)+0.415*range(xlims),min(ylims)+0.275*range(ylims),'draft','fontsize',16,'fontweight','bold');
%iceberg colors
for k = 1:length(highmelt_cmap)
    plot([min(xlims)+0.20*range(xlims) min(xlims)+0.25*range(xlims)],...
        [min(ylims)+0.245*range(ylims)-k*((0.20*range(ylims))/length(highmelt_cmap)) min(ylims)+0.245*range(ylims)-k*((0.20*range(ylims))/length(highmelt_cmap))],...
        '-','linewidth',2*(range(ylims)/range(xlims)),'color',highmelt_cmap((length(highmelt_cmap)+1)-k,:));
end
text(min(xlims)+0.26*range(xlims),min(ylims)+0.245*range(ylims),'50 m yr^{-1}','fontsize',16);
text(min(xlims)+0.26*range(xlims),min(ylims)+0.245*range(ylims)-(length(highmelt_cmap)*4/5)*((0.20*range(ylims))/length(highmelt_cmap)),'10 m yr^{-1}','fontsize',16);
text(min(xlims)+0.26*range(xlims),min(ylims)+0.245*range(ylims)-(0.98*length(highmelt_cmap)*((0.20*range(ylims))/length(highmelt_cmap))),'<1 m yr^{-1}','fontsize',16);
text(min(xlims)+0.23*range(xlims),min(ylims)+0.305*range(ylims),'iceberg','fontsize',16,'fontweight','bold');
text(min(xlims)+0.22*range(xlims),min(ylims)+0.275*range(ylims),'melt rate','fontsize',16,'fontweight','bold');
%add color legend for ocean data
for k = 1:length(Temp_cmap)
    plot([min(xlims)+0.04*range(xlims) min(xlims)+0.09*range(xlims)],...
        [min(ylims)+0.245*range(ylims)-k*((0.20*range(ylims))/length(Temp_cmap)) min(ylims)+0.245*range(ylims)-k*((0.20*range(ylims))/length(Temp_cmap))],...
        '-','color',Temp_cmap((length(Temp_cmap)+1)-k,:)); hold on;
end
text(min(xlims)+0.10*range(xlims),min(ylims)+0.245*range(ylims),[num2str(cmap_add),char(176),'C'],'fontsize',16); 
text(min(xlims)+0.10*range(xlims),min(ylims)+0.245*range(ylims)-(length(Temp_cmap)/2)*((0.20*range(ylims))/length(Temp_cmap)),['0',char(176),'C'],'fontsize',16); 
text(min(xlims)+0.10*range(xlims),min(ylims)+0.245*range(ylims)-(length(Temp_cmap))*((0.20*range(ylims))/length(Temp_cmap)),['-',num2str(cmap_add),char(176),'C'],'fontsize',16);
text(min(xlims)+0.05*range(xlims),min(ylims)+0.305*range(ylims),'ocean','fontsize',16,'fontweight','bold');
text(min(xlims)+0.055*range(xlims),min(ylims)+0.275*range(ylims),'temp.','fontsize',16,'fontweight','bold');

%save figure
[sorted,inds] = sort(leg_ref); mp_sort = mp(inds);
colormap(gca,im_cmap);%make sure the image colormap didn't get accidentally altered
legmap = legend(mp_sort,[char(disp_names(inds))]); set(legmap,'location','northoutside','fontsize',14,'NumColumns',5); 
legmappos = get(legmap,'position'); set(legmap,'position',[0.05 legmappos(2)+0.05 legmappos(3) legmappos(4)]);
gcapos = get(gca,'position'); set(gca,'position',[gcapos(1) 0.09 gcapos(3) gcapos(4)]);
saveas(Tm_mapplot,[figure_path,'Antarctic-iceberg-oceandata-map.eps'],'epsc'); saveas(Tm_mapplot,[figure_path,'Antarctic-iceberg-oceandata-map.png'],'png');


%label the scatterplot & save
figure(Tm_scatterplot); subplot(TF_sub);
set(gca,'fontsize',16); grid on;
set(gca,'xlim',[-0.25 2.25],'ylim',[0 70]); text(2.1,67.5,'a)','fontsize',16);
set(TF_sub,'position',[0.11 0.11 0.40 0.8150]);
%uncomment next 4 lines if you use the plot function to specify symbol colors as a function of draft
% for k = 1:length(depth_cmap)
%     plot([0.1 0.25],[47.5-(k/100) 47.5-(k/100)],'-','color',depth_cmap(k,:)); hold on;
% end
% text(0.275,47,'0 m','fontsize',16); text(0.275,45,'200 m','fontsize',16); text(0.275,39.5,'750 m','fontsize',16);
%next 9 lines should be used if scatterplot function specifies symbol colors as a function of region & symbol size as a function of draft
ylims = get(gca,'ylim'); set(gca,'ylim',[0 max(ylims)]);
rectangle('position',[0.45 max(ylims) - 0.05*((depth_cutoff-50)/150*1.15)*(range(ylims)) 0.525 0.05*((depth_cutoff-50)/150*1.075)*(range(ylims))],'facecolor','w','edgecolor','k');
for j = 1:1:(depth_cutoff-50)/150
    draft_size(j) = round((50+((j-1)*150))*draft_mult)+draft_add;
    yloc(j) = max(ylims) - 0.05*j*(range(ylims)) - 0.005*(range(ylims));
    text(0.6,yloc(j),[num2str((50+((j-1)*150))),' m'],'fontsize',16);
end
scatter(repmat(0.525,size(yloc)),yloc,draft_size,'w','filled','s','markeredgecolor','k'); hold on;
sp_leg = legend(sp,'EAP','EAIS','WAIS','WAP'); set(sp_leg,'location','northwest');
xlabel(['Observed thermal forcing (',char(176),'C)'],'fontsize',16); ylabel('Melt rate (m/yr)','fontsize',16);
% saveas(Tm_scatterplot,[figure_path,'Antarctic-iceberg-meltrate-temp-depth-scatterplots.eps'],'epsc'); saveas(Tm_scatterplot,[figure_path,'Antarctic-iceberg-meltrate-temp-depth-scatterplots.png'],'png');
clear ylims yloc draft_size;
%label and save the SOSE thermal forcing scatterplot
% figure(TmSOSE_scatterplot); 
subplot(TFSOSE_sub);
set(gca,'fontsize',16,'box','on'); grid on;
set(gca,'xlim',[-0.25 2.25],'ylim',[0 70],'yticklabel',[]); text(2.1,67.5,'b)','fontsize',16);
set(TFSOSE_sub,'position',[0.54 0.11 0.40 0.8150]);
% plot ellipses to bound the data
t = linspace(0, 360,1001);
%steeper relationship
xAmplitude = 0.4; yAmplitude = 32; xCenter = 0.15; yCenter = yAmplitude; rotationAngle = -0.6; 
xOriginal = xAmplitude * sind(t) + xCenter; yOriginal = yAmplitude * cosd(t) + yCenter;
transformMatrix = [cosd(rotationAngle), sind(rotationAngle);...
  -sind(rotationAngle), cosd(rotationAngle)];
xyAligned = [xOriginal; yOriginal]'; xyRotated = xyAligned * transformMatrix;
xRotated = xyRotated(:, 1); yRotated = xyRotated(:, 2);
%calculate end points of a line fit to the long axis
xyline = [xOriginal(1),yOriginal(1);xOriginal(501),yOriginal(501)];
xylineRotated = xyline * transformMatrix;
steep_slope = (xylineRotated(1,2) - xylineRotated(2,2))./(xylineRotated(1,1) - xylineRotated(2,1));
fprintf('Approx. increase in melt rate with thermal forcing: %i m/yr per degree C \n',round(steep_slope));
plot(xRotated, yRotated, 'k-.', 'LineWidth', 2,'color',[0.5 0.5 0.5]);
%gradual relationship
xAmplitude = 0.65; yAmplitude = 16; xCenter = 0.70; yCenter = yAmplitude; rotationAngle = -2.25; 
xOriginal = xAmplitude * sind(t) + xCenter; yOriginal = yAmplitude * cosd(t) + yCenter;
transformMatrix = [cosd(rotationAngle), sind(rotationAngle);...
  -sind(rotationAngle), cosd(rotationAngle)];
xyAligned = [xOriginal; yOriginal]'; xyRotated = xyAligned * transformMatrix;
xRotated = xyRotated(:, 1); yRotated = xyRotated(:, 2);
%calculate end points of a line fit to the long axis
xyline = [xOriginal(1),yOriginal(1);xOriginal(501),yOriginal(501)];
xylineRotated = xyline * transformMatrix;
shallow_slope = (xylineRotated(1,2) - xylineRotated(2,2))./(xylineRotated(1,1) - xylineRotated(2,1));
fprintf('Approx. increase in melt rate with thermal forcing: %i m/yr per degree C \n',round(shallow_slope));
plot(xRotated, yRotated, 'k--', 'LineWidth', 2,'color',[0.5 0.5 0.5]);
% ylims = get(gca,'ylim'); set(gca,'ylim',[0 max(ylims)]);
% rectangle('position',[0.175 max(ylims) - 0.05*((depth_cutoff-50)/150*1.15)*(range(ylims)) 0.3 0.05*((depth_cutoff-50)/150*1.05)*(range(ylims))],'facecolor','w','edgecolor','k');
% for j = 1:1:(depth_cutoff-50)/150
%     draft_size(j) = round((50+((j-1)*150))*draft_mult)+draft_add;
%     yloc(j) = max(ylims) - 0.05*j*(range(ylims)) - 0.01*(range(ylims));
%     text(0.275,yloc(j),[num2str((50+((j-1)*150))),' m'],'fontsize',16);
% end
% scatter(repmat(0.225,size(yloc)),yloc,draft_size,'w','filled','s','markeredgecolor','k'); hold on;
% spS_leg = legend(spS,'EAP','EAIS','WAIS','WAP'); set(spS_leg,'location','northwest');
xlabel(['Modeled thermal forcing (',char(176),'C)'],'fontsize',16); %ylabel('Melt rate (m/yr)','fontsize',16);
% saveas(TmSOSE_scatterplot,[figure_path,'Antarctic-iceberg-meltrate-SOSEtemp-depth-scatterplots.eps'],'epsc'); saveas(TmSOSE_scatterplot,[figure_path,'Antarctic-iceberg-meltrate-SOSEtemp-depth-scatterplots.png'],'png');
clear ylims yloc draft_size;
saveas(Tm_scatterplot,[figure_path,'Antarctic-iceberg-meltrate-temp-depth-scatterplots.eps'],'epsc'); saveas(Tm_scatterplot,[figure_path,'Antarctic-iceberg-meltrate-temp-depth-scatterplots.png'],'png');

%add legends to regional maps
%AP
figure(Tm_mapplot_AP); xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
%add color & size legends for iceberg data
scale_yshift = 0.12; scale_ystretch = 0.09; 
% draft_xshift = 0.575; melt_xshift = 0.42; temp_xshift = 0.19; xspan = 0.55; %use if including meltrate color-coding
draft_xshift = 0.35; melt_xshift = 0.20; temp_xshift = melt_xshift; xspan = 0.32; %use if omiting meltrate color-coding
rectangle('position',[max(xlims)-(draft_xshift+0.025)*range(xlims) max(ylims)-(scale_yshift+1.25*scale_ystretch)*range(ylims) xspan*range(xlims) 2.1*scale_ystretch*range(ylims)],'facecolor','w','edgecolor','k'); 
%sizes
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-scale_yshift*range(ylims),round(100*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims),'100 m','fontsize',16); %scaling y-offset = 0.16*range(ylims)
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift+scale_ystretch/2)*range(ylims),round(300*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift+scale_ystretch/2)*range(ylims),'300 m','fontsize',16); %scaling y-offset = 0.11*range(ylims)
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift+scale_ystretch)*range(ylims),round(500*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift+scale_ystretch)*range(ylims),'500 m','fontsize',16); %scaling y-offset = 0.05*range(ylims)
text(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'iceberg','fontsize',16,'fontweight','bold');
text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'draft','fontsize',16,'fontweight','bold');
% %colors
% for k = 1:length(highmelt_cmap)
%     plot([max(xlims)-melt_xshift*range(xlims) max(xlims)-(melt_xshift-0.05)*range(xlims)],...
%         [max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(highmelt_cmap)) max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(highmelt_cmap))],...
%         '-','linewidth',2*(range(ylims)/range(xlims)),'color',highmelt_cmap((length(highmelt_cmap)+1)-k,:));
% end
% text(max(xlims)-(melt_xshift-0.065)*range(xlims),max(ylims)-scale_yshift*range(ylims),'50 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-0.065)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(highmelt_cmap)*4/5)*((scale_ystretch*range(ylims))/length(highmelt_cmap)),'10 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-0.065)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(0.98*length(highmelt_cmap)*((scale_ystretch*range(ylims))/length(highmelt_cmap))),'<1 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-0.045)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'iceberg','fontsize',16,'fontweight','bold');
% text(max(xlims)-(melt_xshift-0.040)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'melt rate','fontsize',16,'fontweight','bold');
%add color legend for ocean temperature
for k = 1:length(Temp_cmap)
    plot([max(xlims)-temp_xshift*range(xlims) max(xlims)-(temp_xshift-0.05)*range(xlims)],...
        [max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(Temp_cmap)) max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(Temp_cmap))],...
        '-','color',Temp_cmap((length(Temp_cmap)+1)-k,:)); hold on;
end
text(max(xlims)-(temp_xshift-0.065)*range(xlims),max(ylims)-scale_yshift*range(ylims),[num2str(cmap_add),char(176),'C'],'fontsize',16); 
text(max(xlims)-(temp_xshift-0.065)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(Temp_cmap)/2)*((scale_ystretch*range(ylims))/length(Temp_cmap)),['0',char(176),'C'],'fontsize',16); 
text(max(xlims)-(temp_xshift-0.065)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(Temp_cmap))*((scale_ystretch*range(ylims))/length(Temp_cmap)),['-',num2str(cmap_add),char(176),'C'],'fontsize',16);
text(max(xlims)-(temp_xshift-0.030)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'ocean','fontsize',16,'fontweight','bold');
text(max(xlims)-(temp_xshift-0.035)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'temp.','fontsize',16,'fontweight','bold');
%THWAITES
figure(Tm_mapplot_TG); xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
%add color & size legends for iceberg data
scale_yshift = 0.87; scale_ystretch = 0.09; 
% draft_xshift = 0.925;  melt_xshift = 0.77; temp_xshift = 0.54; xspan = 0.55; %use if including meltrate color-coding
draft_xshift = 0.925; melt_xshift = 0.775; temp_xshift = melt_xshift; xspan = 0.32; %use if omiting meltrate color-coding
rectangle('position',[max(xlims)-(draft_xshift+0.025)*range(xlims) max(ylims)-(scale_yshift+1.25*scale_ystretch)*range(ylims) xspan*range(xlims) 2.1*scale_ystretch*range(ylims)],'facecolor','w','edgecolor','k'); 
%sizes
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-scale_yshift*range(ylims),round(100*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims),'100 m','fontsize',16); %scaling y-offset = 0.16*range(ylims)
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift+scale_ystretch/2)*range(ylims),round(300*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift+scale_ystretch/2)*range(ylims),'300 m','fontsize',16); %scaling y-offset = 0.11*range(ylims)
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift+scale_ystretch)*range(ylims),round(500*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift+scale_ystretch)*range(ylims),'500 m','fontsize',16); %scaling y-offset = 0.05*range(ylims)
text(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'iceberg','fontsize',16,'fontweight','bold');
text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'draft','fontsize',16,'fontweight','bold');
% %colors
% for k = 1:length(highmelt_cmap)
%     plot([max(xlims)-melt_xshift*range(xlims) max(xlims)-(melt_xshift-0.05)*range(xlims)],...
%         [max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(highmelt_cmap)) max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(highmelt_cmap))],...
%         '-','linewidth',2*(range(ylims)/range(xlims)),'color',highmelt_cmap((length(highmelt_cmap)+1)-k,:));
% end
% text(max(xlims)-(melt_xshift-0.065)*range(xlims),max(ylims)-scale_yshift*range(ylims),'50 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-0.065)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(highmelt_cmap)*4/5)*((scale_ystretch*range(ylims))/length(highmelt_cmap)),'10 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-0.065)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(0.98*length(highmelt_cmap)*((scale_ystretch*range(ylims))/length(highmelt_cmap))),'<1 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-0.045)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'iceberg','fontsize',16,'fontweight','bold');
% text(max(xlims)-(melt_xshift-0.040)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'melt rate','fontsize',16,'fontweight','bold');
%add color legend for ocean temperature
for k = 1:length(Temp_cmap)
    plot([max(xlims)-temp_xshift*range(xlims) max(xlims)-(temp_xshift-0.05)*range(xlims)],...
        [max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(Temp_cmap)) max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(Temp_cmap))],...
        '-','color',Temp_cmap((length(Temp_cmap)+1)-k,:)); hold on;
end
text(max(xlims)-(temp_xshift-0.065)*range(xlims),max(ylims)-scale_yshift*range(ylims),[num2str(cmap_add),char(176),'C'],'fontsize',16); 
text(max(xlims)-(temp_xshift-0.065)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(Temp_cmap)/2)*((scale_ystretch*range(ylims))/length(Temp_cmap)),['0',char(176),'C'],'fontsize',16); 
text(max(xlims)-(temp_xshift-0.065)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(Temp_cmap))*((scale_ystretch*range(ylims))/length(Temp_cmap)),['-',num2str(cmap_add),char(176),'C'],'fontsize',16);
text(max(xlims)-(temp_xshift-0.030)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'ocean','fontsize',16,'fontweight','bold');
text(max(xlims)-(temp_xshift-0.035)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'temp.','fontsize',16,'fontweight','bold');
%MERTZ
figure(Tm_mapplot_MI); xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
%add color & size legends for iceberg data
scale_yshift = 0.75; scale_ystretch = 0.18; 
% draft_xshift = 0.925; melt_xshift = 0.81; temp_xshift = 0.62; xspan = 0.45;%use if including meltrate color-coding
draft_xshift = 0.95; melt_xshift = 0.815; temp_xshift = melt_xshift; xspan = 0.27; %use if omiting meltrate color-coding
rectangle('position',[max(xlims)-(draft_xshift+0.025)*range(xlims) max(ylims)-(scale_yshift+1.25*scale_ystretch)*range(ylims) xspan*range(xlims) 2.1*scale_ystretch*range(ylims)],'facecolor','w','edgecolor','k'); 
%sizes
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-scale_yshift*range(ylims),round(100*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims),'100 m','fontsize',16); %scaling y-offset = 0.16*range(ylims)
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift+scale_ystretch/2)*range(ylims),round(300*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift+scale_ystretch/2)*range(ylims),'300 m','fontsize',16); %scaling y-offset = 0.11*range(ylims)
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift+scale_ystretch)*range(ylims),round(500*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift+scale_ystretch)*range(ylims),'500 m','fontsize',16); %scaling y-offset = 0.05*range(ylims)
text(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'iceberg','fontsize',16,'fontweight','bold');
text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'draft','fontsize',16,'fontweight','bold');
% %colors
% for k = 1:length(highmelt_cmap)
%     plot([max(xlims)-melt_xshift*range(xlims) max(xlims)-(melt_xshift-(xspan/10))*range(xlims)],...
%         [max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(highmelt_cmap)) max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(highmelt_cmap))],...
%         '-','linewidth',2*(range(ylims)/range(xlims)),'color',highmelt_cmap((length(highmelt_cmap)+1)-k,:));
% end
% text(max(xlims)-(melt_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims),'50 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(highmelt_cmap)*4/5)*((scale_ystretch*range(ylims))/length(highmelt_cmap)),'10 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(0.98*length(highmelt_cmap)*((scale_ystretch*range(ylims))/length(highmelt_cmap))),'<1 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-(xspan/10)+0.015)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'iceberg','fontsize',16,'fontweight','bold');
% text(max(xlims)-(melt_xshift-(xspan/10)+0.020)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'melt rate','fontsize',16,'fontweight','bold');
%add color legend for ocean temperature
for k = 1:length(Temp_cmap)
    plot([max(xlims)-temp_xshift*range(xlims) max(xlims)-(temp_xshift-(xspan/10))*range(xlims)],...
        [max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(Temp_cmap)) max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(Temp_cmap))],...
        '-','color',Temp_cmap((length(Temp_cmap)+1)-k,:)); hold on;
end
text(max(xlims)-(temp_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims),[num2str(cmap_add),char(176),'C'],'fontsize',16); 
text(max(xlims)-(temp_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(Temp_cmap)/2)*((scale_ystretch*range(ylims))/length(Temp_cmap)),['0',char(176),'C'],'fontsize',16); 
text(max(xlims)-(temp_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(Temp_cmap))*((scale_ystretch*range(ylims))/length(Temp_cmap)),['-',num2str(cmap_add),char(176),'C'],'fontsize',16);
text(max(xlims)-(temp_xshift-(xspan/10)+0.025)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'ocean','fontsize',16,'fontweight','bold');
text(max(xlims)-(temp_xshift-(xspan/10)+0.02)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'temp.','fontsize',16,'fontweight','bold');
%TOTTEN
figure(Tm_mapplot_TI); xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
%add color & size legends for iceberg data
scale_yshift = 0.75; scale_ystretch = 0.18; 
% draft_xshift = 0.475; melt_xshift = 0.36; temp_xshift = 0.17; xspan = 0.45;%use if including meltrate color-coding
draft_xshift = 0.275; melt_xshift = 0.155; temp_xshift = melt_xshift; xspan = 0.25; %use if omiting meltrate color-coding
rectangle('position',[max(xlims)-(draft_xshift+0.025)*range(xlims) max(ylims)-(scale_yshift+1.25*scale_ystretch)*range(ylims) xspan*range(xlims) 2.1*scale_ystretch*range(ylims)],'facecolor','w','edgecolor','k'); 
%sizes
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-scale_yshift*range(ylims),round(100*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims),'100 m','fontsize',16); %scaling y-offset = 0.16*range(ylims)
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift+scale_ystretch/2)*range(ylims),round(300*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift+scale_ystretch/2)*range(ylims),'300 m','fontsize',16); %scaling y-offset = 0.11*range(ylims)
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift+scale_ystretch)*range(ylims),round(500*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift+scale_ystretch)*range(ylims),'500 m','fontsize',16); %scaling y-offset = 0.05*range(ylims)
text(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'iceberg','fontsize',16,'fontweight','bold');
text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'draft','fontsize',16,'fontweight','bold');
% %colors
% for k = 1:length(highmelt_cmap)
%     plot([max(xlims)-melt_xshift*range(xlims) max(xlims)-(melt_xshift-(xspan/10))*range(xlims)],...
%         [max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(highmelt_cmap)) max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(highmelt_cmap))],...
%         '-','linewidth',2*(range(ylims)/range(xlims)),'color',highmelt_cmap((length(highmelt_cmap)+1)-k,:));
% end
% text(max(xlims)-(melt_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims),'50 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(highmelt_cmap)*4/5)*((scale_ystretch*range(ylims))/length(highmelt_cmap)),'10 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(0.98*length(highmelt_cmap)*((scale_ystretch*range(ylims))/length(highmelt_cmap))),'<1 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-(xspan/10)+0.015)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'iceberg','fontsize',16,'fontweight','bold');
% text(max(xlims)-(melt_xshift-(xspan/10)+0.020)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'melt rate','fontsize',16,'fontweight','bold');
%add color legend for ocean temperature
for k = 1:length(Temp_cmap)
    plot([max(xlims)-temp_xshift*range(xlims) max(xlims)-(temp_xshift-(xspan/10))*range(xlims)],...
        [max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(Temp_cmap)) max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(Temp_cmap))],...
        '-','color',Temp_cmap((length(Temp_cmap)+1)-k,:)); hold on;
end
text(max(xlims)-(temp_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims),[num2str(cmap_add),char(176),'C'],'fontsize',16); 
text(max(xlims)-(temp_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(Temp_cmap)/2)*((scale_ystretch*range(ylims))/length(Temp_cmap)),['0',char(176),'C'],'fontsize',16); 
text(max(xlims)-(temp_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(Temp_cmap))*((scale_ystretch*range(ylims))/length(Temp_cmap)),['-',num2str(cmap_add),char(176),'C'],'fontsize',16);
text(max(xlims)-(temp_xshift-(xspan/10)+0.025)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'ocean','fontsize',16,'fontweight','bold');
text(max(xlims)-(temp_xshift-(xspan/10)+0.02)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'temp.','fontsize',16,'fontweight','bold');
%FILCHNER
figure(Tm_mapplot_FI); xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
%add color & size legends for iceberg data
scale_yshift = 0.87; scale_ystretch = 0.09; 
% draft_xshift = 0.55; melt_xshift = 0.40; temp_xshift = 0.18; xspan = 0.54;%use if including meltrate color-coding
draft_xshift = 0.33; melt_xshift = 0.18; temp_xshift = melt_xshift; xspan = 0.30; %use if omiting meltrate color-coding
rectangle('position',[max(xlims)-(draft_xshift+0.025)*range(xlims) max(ylims)-(scale_yshift+1.25*scale_ystretch)*range(ylims) xspan*range(xlims) 2.1*scale_ystretch*range(ylims)],'facecolor','w','edgecolor','k'); 
%sizes
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-scale_yshift*range(ylims),round(100*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims),'100 m','fontsize',16); %scaling y-offset = 0.16*range(ylims)
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift+scale_ystretch/2)*range(ylims),round(300*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift+scale_ystretch/2)*range(ylims),'300 m','fontsize',16); %scaling y-offset = 0.11*range(ylims)
scatter(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift+scale_ystretch)*range(ylims),round(500*draft_mult+draft_add),'w','filled','s','markeredgecolor','k'); text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift+scale_ystretch)*range(ylims),'500 m','fontsize',16); %scaling y-offset = 0.05*range(ylims)
text(max(xlims)-draft_xshift*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'iceberg','fontsize',16,'fontweight','bold');
text(max(xlims)-(draft_xshift-0.015)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'draft','fontsize',16,'fontweight','bold');
% %colors
% for k = 1:length(highmelt_cmap)
%     plot([max(xlims)-melt_xshift*range(xlims) max(xlims)-(melt_xshift-(xspan/10))*range(xlims)],...
%         [max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(highmelt_cmap)) max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(highmelt_cmap))],...
%         '-','linewidth',2*(range(ylims)/range(xlims)),'color',highmelt_cmap((length(highmelt_cmap)+1)-k,:));
% end
% text(max(xlims)-(melt_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims),'50 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(highmelt_cmap)*4/5)*((scale_ystretch*range(ylims))/length(highmelt_cmap)),'10 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(0.98*length(highmelt_cmap)*((scale_ystretch*range(ylims))/length(highmelt_cmap))),'<1 m yr^{-1}','fontsize',16);
% text(max(xlims)-(melt_xshift-(xspan/10)+0.015)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'iceberg','fontsize',16,'fontweight','bold');
% text(max(xlims)-(melt_xshift-(xspan/10)+0.020)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'melt rate','fontsize',16,'fontweight','bold');
%add color legend for ocean temperature
for k = 1:length(Temp_cmap)
    plot([max(xlims)-temp_xshift*range(xlims) max(xlims)-(temp_xshift-(xspan/10))*range(xlims)],...
        [max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(Temp_cmap)) max(ylims)-scale_yshift*range(ylims)-k*((scale_ystretch*range(ylims))/length(Temp_cmap))],...
        '-','color',Temp_cmap((length(Temp_cmap)+1)-k,:)); hold on;
end
text(max(xlims)-(temp_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims),[num2str(cmap_add),char(176),'C'],'fontsize',16); 
text(max(xlims)-(temp_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(Temp_cmap)/2)*((scale_ystretch*range(ylims))/length(Temp_cmap)),['0',char(176),'C'],'fontsize',16); 
text(max(xlims)-(temp_xshift-(xspan/10)-0.015)*range(xlims),max(ylims)-scale_yshift*range(ylims)-(length(Temp_cmap))*((scale_ystretch*range(ylims))/length(Temp_cmap)),['-',num2str(cmap_add),char(176),'C'],'fontsize',16);
text(max(xlims)-(temp_xshift-(xspan/10)+0.025)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/1.5)*range(ylims),'ocean','fontsize',16,'fontweight','bold');
text(max(xlims)-(temp_xshift-(xspan/10)+0.02)*range(xlims),max(ylims)-(scale_yshift-scale_ystretch/3)*range(ylims),'temp.','fontsize',16,'fontweight','bold');

%save regional maps
saveas(Tm_mapplot_AP,[figure_path,'AP_iceberg-oceandata-map.eps'],'epsc'); saveas(Tm_mapplot_AP,[figure_path,'AP_iceberg-oceandata-map.png'],'png');
saveas(Tm_mapplot_TG,[figure_path,'TG_iceberg-oceandata-map.eps'],'epsc'); saveas(Tm_mapplot_TG,[figure_path,'TG_iceberg-oceandata-map.png'],'png');
saveas(Tm_mapplot_MI,[figure_path,'MI_iceberg-oceandata-map.eps'],'epsc'); saveas(Tm_mapplot_MI,[figure_path,'MI_iceberg-oceandata-map.png'],'png');
saveas(Tm_mapplot_TI,[figure_path,'TI_iceberg-oceandata-map.eps'],'epsc'); saveas(Tm_mapplot_TI,[figure_path,'TI_iceberg-oceandata-map.png'],'png');
saveas(Tm_mapplot_FI,[figure_path,'FI_iceberg-oceandata-map.eps'],'epsc'); saveas(Tm_mapplot_FI,[figure_path,'FI_iceberg-oceandata-map.png'],'png');

disp('plots completed!');

%% export summary data on iceberg melt rates & ocean observations to tables

%create dummy matrices to hold the summary data
disp_name = {}; 
dateo = []; datef = [];
noobs = [];
xmin = []; xmax = []; ymin = []; ymax = [];
d = []; dmin = []; dmax = [];
m = []; mmin = []; mmax = [];
ind = 1;

%loop through the iceberg data, showing info, & exporting to a csv
disp('printing out stats');
for j = 1:length(melt)
%     figure; title(string(melt(j).dispname));
    %find unique date ranges
    date_pairs = [melt(j).to,melt(j).tf];
    [unique_dates,unique_inds,date_refs] = unique(date_pairs,'rows');
    
    %create filters for melt rate vs draft stats
    %identify outliers as anomalously high values
    outliers = isoutlier(melt(j).m); 
    m_filtered = 365*melt(j).m; 
    filter_ref = find(outliers == 1 & m_filtered > nanmedian(m_filtered)+6*1.4826*mad(m_filtered,1));
    m_filtered(filter_ref) = NaN;
    %create geographic filter (for Thwaites & Mertz)
    region = ones(size(melt(j).m));
    if strmatch('Thwaites',string(melt(j).dispname))
        East_daterefs = [1,3,5]; %East_dates = unique_dates(East_daterefs,:)
        for k = 1:length(East_daterefs); region(date_refs == East_daterefs(k)) = 0; end
        clear East_date*;
    elseif strmatch('Mertz',string(melt(j).dispname))
        East_daterefs = 4; %East_dates = unique_dates(East_daterefs,:)
        for k = 1:length(East_daterefs); region(date_refs == East_daterefs(k)) = 0; end
        clear East_date*;
    end
%     %plot distributions of meltrate vs draft 
%     histogram(m_filtered./melt(j).d,'FaceColor','k'); hold on;
%     histogram(m_filtered(region==1)./melt(j).d(region==1),'FaceColor','w'); hold on;
    
    %loop through dates, summarizing data
    for k = 1:max(date_refs)
        to(k,:) = convert_from_decimaldate(melt(j).to(unique_inds(k)));
        tf(k,:) = convert_from_decimaldate(melt(j).tf(unique_inds(k)));
        x_range(k,:) = [round(min(melt(j).x(date_refs == k))),round(max(melt(j).x(date_refs == k)))];
        y_range(k,:) = [round(min(melt(j).y(date_refs == k))),round(max(melt(j).y(date_refs == k)))];
        obs(k,1) = length(find(date_refs == k));
        draft(k,1) = nanmedian(melt(j).d(date_refs == k));
        draft_range(k,:) = [min(melt(j).d(date_refs == k)),max(melt(j).d(date_refs == k))];
        meltrate(k,1) = 365*nanmedian(melt(j).m(date_refs == k));
        meltrate_range(k,:) = 365*[min(melt(j).m(date_refs == k)),max(melt(j).m(date_refs == k))];
        %print summary stats to the command window
        fprintf('%s, %s, %s, %i, (%7.0f:%7.0f), (%7.0f:%7.0f), %4.1f (%4.1f:%4.1f), %3.1f (%3.2f:%3.1f) & melt/draft = %4.3f (site median,mad = %4.3f, %4.3f) \n',...
            string(melt(j).dispname),... %site name
            [to(k,1:4),'/',to(k,5:6),'/',to(k,7:8)],[tf(k,1:4),'/',tf(k,5:6),'/',tf(k,7:8)], obs(k),... %dates
            max(x_range(k,:)),min(x_range(k,:)),max(y_range(k,:)),min(y_range(k,:)), draft(k),...
            min(draft_range(k,:)),max(draft_range(k,:)),meltrate(k),min(meltrate_range(k,:)),max(meltrate_range(k,:)),...
            nanmedian(m_filtered(date_refs == k)./melt(j).d(date_refs == k)),...
            nanmedian(m_filtered(region==1)./melt(j).d(region==1)),mad(m_filtered(region==1)./melt(j).d(region==1),1));
%         histogram(m_filtered(date_refs == k)./melt(j).d(date_refs == k)); hold on;
        %add summary stats to variables that will be exported to a csv
        disp_name(ind,1) = melt(j).dispname;
        dateo = [dateo; to(k,1:4),'/',to(k,5:6),'/',to(k,7:8)]; datef = [datef; tf(k,1:4),'/',tf(k,5:6),'/',tf(k,7:8)];
        noobs(ind,1) = obs(k,1);
        xmin(ind,1) = x_range(k,1); xmax(ind,1) = x_range(k,2); ymin(ind,1) = y_range(k,1); ymax(ind,1) = y_range(k,2);
        d(ind,1) = round(draft(k,1)); dmin(ind,1) = round(draft_range(k,1)); dmax(ind,1) = round(draft_range(k,2));
        m(ind,1) = round(meltrate(k,1),1); mmin(ind,1) = round(meltrate_range(k,1),1); mmax(ind,1) = round(meltrate_range(k,2),1);
        
        %advance to the next date pair
        ind = ind+1;
        clear to tf *_range obs draft* meltrate*;
    end
    
    %advance to the next site
    clear date_pairs unique_* date_refs outliers filter_ref m_filtered;
end
%save iceberg summary data to a table
T=table(disp_name, dateo, datef, noobs, xmin, xmax, ymin, ymax, d, dmin, dmax, m, mmin, mmax);
column_names = ["Site Name", "Early Date", "Late Date", "Iceberg Count",...
    "X minimum (PS m)", "X maximum (PS m)", "Y minimum (PS m)", "Y maximum (PS m)",...
    "Draft (m b.s.l.)", "Draft minimum (m b.s.l.)", "Draft maximum (m b.s.l.)",...
    "Meltrate (m/yr)", "Meltrate minimum (m/yr)", "Meltrate maximum (m/yr)"];
T.Properties.VariableNames = column_names;
writetable(T,[figure_path,'/','Antarctic_iceberg_summary-stats.csv']);
clear T;

%set-up dummy matrices for ocean temp data
clear disp_name ind; disp_name = {}; 
ind = 1;

%loop through the ocean data, showing info, & exporting to a csv
for j = 1:length(melt)
    if size(melt(j).oceant,1) > 0
        disp_name(ind,1) = melt(j).dispname;
        obs_all(ind,1) = size(melt(j).oceant,1);
        obs_dated(ind,1) = size(melt(j).oceantavg,2);
        Tmedian(ind,1) = round(nanmedian(melt(j).oceanTavg),1); 
        Tmin(ind,1) = round(min(melt(j).oceanTavg),1); Tmax(ind,1) = round(max(melt(j).oceanTavg),1);
        TFmedian(ind,:) = round(nanmedian(melt(j).oceanTavg-melt(j).oceanTfp),1); 
        TFmin(ind,1) = round(min(melt(j).oceanTavg-melt(j).oceanTfp),1); TFmax(ind,1) = round(max(melt(j).oceanTavg-melt(j).oceanTfp),1);
        ind = ind + 1;
    end
end
%save iceberg summary data to a table
T=table(disp_name, obs_all, obs_dated, Tmedian, Tmin, Tmax, TFmedian, TFmin, TFmax);
column_names = ["Site Name", "Profile Count", "Unique Date Count",...
    "Temperature (Celsius)", "Temperature minimum (Celsius)", "Temperature maximum (Celsius)",...
    "Thermal Forcing (Celsius)", "Thermal Forcing minimum (Celsius)", "Thermal Forcing maximum (Celsius)"];
T.Properties.VariableNames = column_names;
writetable(T,[figure_path,'/','Antarctic_oceantemp_summary-stats.csv']);
clear T;

