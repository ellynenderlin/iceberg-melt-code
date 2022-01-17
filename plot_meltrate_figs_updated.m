clearvars; close all; drawnow;
addpath('/users/ellynenderlin/mfiles','/users/ellynenderlin/mfiles/general','/users/ellynenderlin/mfiles/general/cmocean');

cd /Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/
plot_yrs = []; avgx = []; avgy = []; region = []; warning off;
rho_sw = 1026; %sea water density in kg m^-3

%make colormaps
yrs = [2013:1:2020]; %plot_marker = colormap(parula(length(yrs)+1)); plot_marker = plot_marker(1:length(yrs),:);
% plot_marker = [127,59,8;179,88,6;224,130,20;253,184,99;254,224,182;216,218,235;178,171,210;128,115,172;84,39,136;45,0,75]/255; %orange-purple colormap
plot_marker = cmocean('dense',9); plot_marker = plot_marker(2:end,:); %colormap for emphasizing different years
lowmelt_cmap = flipud(colormap(hot(40))); highmelt_cmap = flipud(colormap(hot(100)));
% plot_cmap = colormap(parula(15)); %colormap for emphasizing different locations
close all;

%if subplots for each region are organized 2 rows by 8 columns
region = [{'Edgeworth-LarsenA'},{'Crane-LarsenB'},{'RonneW'},{'Filchner'},{'Amery'},{'Totten'},{'MertzW'},{'Thwaites'},{'Eltanin'},{'Seller-Bugge'},{'Heim-Marguerite'},{'Widdowson-Biscoe'},{'Cadman-Palmer'},{'Blanchard-Danco'},{'Leonardo-Danco'}];
disp_names = [{'i) Edgeworth'},{'j) Crane'},{'k) Ronne'},{'l) Filchner'},{'m) Amery'},{'n) Totten'},{'o) Mertz'},{'h) Thwaites'},{'g) Ferrigno'},{'f) Seller'},{'e) Heim'},{'d) Widdowson'},{'c) Cadman'},{'b) Blanchard'},{'a) Leonardo'}];
leg_names = [{'Edgeworth'},{'Crane'},{'Ronne'},{'Filchner'},{'Amery'},{'Totten'},{'Mertz'},{'Thwaites'},{'Ferrigno'},{'Seller'},{'Heim'},{'Widdowson'},{'Cadman'},{'Blanchard'},{'Leonardo'}];
% marker = ['d','d','d','d','d','d','d','s','s','s','s','s','s','s','s']; %set-up the different marker styles for E and W
marker = ['s','s','s','s','s','s','s','s','s','s','s','s','s','s','s']; %set-up the different marker styles for E and W
plot_loc = [2,4,6,8,10,12,14,15,13,11,9,7,5,3,1];
%if subplots for each region are organized 8 rows by 2 columns
% region = [{'Leonardo'},{'LarsenA'},{'Danco'},{'LarsenB'},{'Palmer'},{'Ronne'},{'Biscoe'},{'Filchner'},{'Marguerite'},{'Amery'},{'Bugge'},{'Totten'},{'Eltanin'},{'Mertz'},{'Thwaites'}];
% disp_names = [{'Leonardo'},{'Edgeworth'},{'Blanchard'},{'Crane'},{'Cadman'},{'Ronne'},{'Widdowson'},{'Filchner'},{'Heim'},{'Christensen'},{'Seller'},{'Totten'},{'Eltanin'},{'Mertz'},{'Thwaites'}];
% marker = ['h','s','h','s','h','d','h','d','h','o','h','o','p','o','p'];



%set-up a location map
cd /Users/ellynenderlin/Research/miscellaneous/RAMP
[A,S] = geotiffread('Antarctic_RAMP_image_v2_1km.tif');
IM.x = S.XWorldLimits(1)+0.5*S.CellExtentInWorldX:S.CellExtentInWorldX:S.XWorldLimits(2)-0.5*S.CellExtentInWorldX;
IM.y = S.YWorldLimits(2)-0.5*S.CellExtentInWorldY:-S.CellExtentInWorldY:S.YWorldLimits(1)+0.5*S.CellExtentInWorldY;
IM.z=single(A);
clear A S;
%set-up a map
figureA=figure; set(gcf,'position',[450 50 800 800]);
im_cmap = colormap(gray(10001)); im_cmap(1,:) = [1 1 1];
imagesc(IM.x,IM.y,IM.z); colormap(gca,im_cmap); hold on; axis xy equal
%set-up subplots for graphs
figureB = figure; set(gcf,'position',[50 400 1600 600]);
sub1b = subplot(1,3,1); sub2b = subplot(1,3,2); sub3b = subplot(1,3,3);
figureC = figure; set(gcf,'position',[50 50 800 1000]);

%create dummy matrices to fill with meltwater flux & submerged area for E
%vs W Antarctica in order to fit trendlines to the 2 datasets
WdVdt = []; WAsub = []; EdVdt = []; EAsub = [];

%set-up the plots for the combined datasets
cd /Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt
for i = 1:size(region,2)
    cd_to_dir = ['cd ',char(region(i))]; eval(cd_to_dir);
    meltinfo = dir('*iceberg_meltinfo.txt'); landsats = dir('LC*PS.TIF');
    date_o = []; xcoord_o = []; ycoord_o = [];
    date_f = []; xcoord_f = []; ycoord_f = [];
    flux = []; sub_area = []; meltrate = []; keeld = [];
    figureE = figure; set(gcf,'position',[450 450 800 300]);
    for j = 1:length(meltinfo)
        plot_yrs = [str2num(meltinfo(j).name(end-49:end-46)) str2num(meltinfo(j).name(end-34:end-31))];
        load_data = ['M=dlmread(''',meltinfo(j).name,''');']; eval(load_data);
        
        %identify data with clear issues
        bad_data = find(M(:,18)<0); M(bad_data,:) = [];
        
        %pull variables
        dt = M(:,1); %time
        xo = M(:,2); yo = M(:,3); zo = M(:,4); rhoo = M(:,5); Vo = M(:,6); %initial locations, median elev, density, volume
        xf = M(:,7); yf = M(:,8); zf = M(:,9); rhof = M(:,10); Vf = M(:,11); %same as above but final
        coregzo = M(:,12); coregzf = M(:,13);
        dz = M(:,14); dz_sigma = M(:,15);
        dVdt = M(:,16); dVdt_uncert = M(:,17);

        %discovered some odd issue with the exported drafts, recalculating those and the submerged areas
        draft = (nanmean([rhoo rhof],2)./(repmat(rho_sw,length(nanmean([rhoo rhof],2)),1)-nanmean([rhoo rhof],2))).*nanmean([zo zf],2); %draft = M(:,18); 
        draft_uncert = M(:,19);
        Asurf = M(:,20); Asurf_uncert = M(:,21);
        lat_area = M(:,22) - Asurf; perim = lat_area./draft; clear lat_area; lat_area = perim.*draft; Asub = lat_area + Asurf; clear lat_area perim; %Asub = M(:,22); 
        Asub_uncert = M(:,23);
        avgx(i) = nanmean(xo); avgy(i) = nanmean(yo); %this determines where symbols get plotted (i.e. study site locations)
        date_o = [date_o; repmat(str2num(meltinfo(j).name(4:11)),size(xo))]; xcoord_o = [xcoord_o; xo]; ycoord_o = [ycoord_o; yo];
        date_f = [date_f; repmat(str2num(meltinfo(j).name(13:20)),size(xf))]; xcoord_f = [xcoord_f; xf]; ycoord_f = [ycoord_f; yf];
        flux = [flux; dVdt]; sub_area = [sub_area; Asub]; meltrate = [meltrate; (dVdt./Asub)]; keeld = [keeld; draft];
        figure(figureA); colormap(gca,'gray');
        plot(avgx(i),avgy(i),[marker(i),'k'],'markerfacecolor','w','markeredgecolor','k','markersize',16); hold on;
%         if j == length(meltinfo)
%             if mod(plot_loc(i),2) == 0 %west Antarctica
%                 text(avgx(i)-1000000,avgy(i),disp_names(i),'fontsize',14,'color','w');
%                 WdVdt = [WdVdt; dVdt]; WAsub = [WAsub; Asub];
%             else %east Antarctica
%                  text(avgx(i)+50000,avgy(i),disp_names(i),'fontsize',14,'color','w');
%                  EdVdt = [EdVdt; dVdt]; EAsub = [EAsub; Asub];
%             end
%         end
        if ~isempty(strmatch(marker(i),'p')) || ~isempty(strmatch(marker(i),'h')); symbol_size = 10; else symbol_size = 8; end
        
        %export avg data for each date all to one table
        cd ..
        table_data = [str2num(meltinfo(j).name(end-49:end-46)) str2num(meltinfo(j).name(end-34:end-31)) nanmedian(dVdt) nanmedian(dVdt_uncert) nanmedian(draft) nanmedian(draft_uncert) nanmedian(Asub) nanmedian(Asub_uncert) nanmedian(dVdt./Asub) nanmedian(abs(dVdt).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2))];  
        if i == 1 && j == 1
            dlmwrite('Antarctic-iceberg-meltinfo.txt',table_data,'delimiter','\t','precision',8);
        else
            dlmwrite('Antarctic-iceberg-meltinfo.txt',table_data,'-append','delimiter','\t','precision',8);
        end
        clear table_data
        cd_to_dir = ['cd ',char(region(i))]; eval(cd_to_dir);
        
        %multi-panel subplots of all data
        figure(figureB);
        subplot(sub1b);
        %color-coded to region
%         errorbar(Asub,dVdt/86400,dVdt_uncert/86400,dVdt_uncert/86400,Asub_uncert,Asub_uncert,[marker(i),'k'],'markerfacecolor',plot_cmap(i,:),'markersize',symbol_size); hold on;
%         pl(i) = plot(Asub,dVdt/86400,[marker(i),'k'],'markerfacecolor',plot_cmap(i,:),'markersize',symbol_size); hold on;
%         subplot(sub2b);
%         errorbar(draft,dVdt./Asub,abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,[marker(i),'k'],'markerfacecolor',plot_cmap(i,:),'markersize',symbol_size); hold on;
%         plot(draft,dVdt./Asub,[marker(i),'k'],'markerfacecolor',plot_cmap(i,:),'markersize',symbol_size); hold on;
%         if i < 10
%             subplot(sub3b);
%             errorbar(draft,dVdt./Asub,abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,[marker(i),'k'],'markerfacecolor',plot_cmap(i,:),'markersize',symbol_size); hold on;
%             plot(draft,dVdt./Asub,[marker(i),'k'],'markerfacecolor',plot_cmap(i,:),'markersize',symbol_size); hold on;
%         end
        %color-coded to year
        errorbar(Asub,dVdt/86400,dVdt_uncert/86400,dVdt_uncert/86400,Asub_uncert,Asub_uncert,[marker(i),'k'],'markerfacecolor',plot_marker(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markersize',symbol_size); hold on;
%         if j == 1
            pl(i) = plot(Asub,dVdt/86400,[marker(i),'k'],'markerfacecolor','w','markersize',symbol_size); hold on;
%         end
%         plot(Asub,365*dVdt,[marker(i),'k'],'markerfacecolor',plot_marker(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markersize',symbol_size); hold on;
        subplot(sub2b);
        errorbar(draft,365*dVdt./Asub,365*abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),365*abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,[marker(i),'k'],'markerfacecolor','w','markersize',symbol_size); hold on;
        plot(draft,365*dVdt./Asub,[marker(i),'k'],'markerfacecolor','w','markersize',symbol_size); hold on;
        if i < 10
            subplot(sub3b);
            errorbar(draft,365*dVdt./Asub,365*abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),365*abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,[marker(i),'k'],'markerfacecolor','w','markersize',symbol_size); hold on;
            plot(draft,365*dVdt./Asub,[marker(i),'k'],'markerfacecolor','w','markersize',symbol_size); hold on;
        end
        
        %subplots for each region but all in one figure window
        figure(figureC);
%         if i>=5; subplot(3,5,i+1); else; subplot(3,5,i); end
        if j == 1
            subpl = subplot(8,2,plot_loc(i));
        else
            subplot(subpl);
        end
        %color-coded according to region w/ edges color-coded by date
%         errorbar(draft,dVdt./Asub,abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,[marker(i),'k'],'markerfacecolor',plot_cmap(i,:),'markeredgecolor',plot_marker(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markersize',symbol_size); hold on;
%         plot(draft,dVdt./Asub,[marker(i),'k'],'markerfacecolor',plot_cmap(i,:),'markeredgecolor',plot_marker(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markersize',symbol_size,'linewidth',1.5); hold on;
        %add a legend to lower right plot
        plotpos = get(gca,'position');
        if plot_loc(i) == 14
            if j == 1
                for k = 1:length(yrs)
                    yrpl(k) = plot(draft(1),365*dVdt(1)./Asub(1),'sk','markerfacecolor',plot_marker(yrs(k)-(min(yrs)-1),:)); hold on;
                end
            end
            set(gca,'xlim',[0 400],'xtick',[0:100:400],'xticklabel',[0:100:400]);
        end
        %plot face color color-coded by date
        errorbar(draft,365*dVdt./Asub,365*abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),365*abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,[marker(i),'k'],'markerfacecolor',plot_marker(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markeredgecolor','k','markersize',symbol_size); hold on;
        plot(draft,365*dVdt./Asub,[marker(i),'k'],'markerfacecolor',plot_marker(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markeredgecolor','k','markersize',symbol_size); hold on;
        if avgx(i) > -1.5e6
            set(gca,'xlim',[0 400],'xtick',[0:100:400],'xticklabel',[],...
                    'ylim',[0 9],'ytick',[0:4:9],'yticklabel',[0:4:9],'fontsize',16); grid on;
        else
            if avgx(i) > -2.45e6 && avgy(i) > 1.2e6
                set(gca,'xlim',[0 400],'xtick',[0:100:400],'xticklabel',[],...
                    'ylim',[0 22],'ytick',[0:10:22],'yticklabel',[0:10:22],'fontsize',16); grid on;
            else
                set(gca,'xlim',[0 800],'xtick',[0:200:800],'xticklabel',[],...
                    'ylim',[0 70],'ytick',[0:30:70],'yticklabel',[0:30:70],'fontsize',16); grid on;
            end
        end
        
        %add axes labels
        if plot_loc(i) == 14
            if j == length(meltinfo)
                legsub = legend(yrpl,num2str(yrs')); set(legsub,'location','southoutside','NumColumns',4,'fontsize',16);
                set(gca,'xlim',[0 400],'xtick',[0:100:400],'xticklabel',[0:100:400]); grid on;
                xlabel('Draft (m b.s.l.)','fontsize',20);
            end
        end
        if plot_loc(i) == 15
            xlabel('Draft (m b.s.l.)','fontsize',20); ylb = ylabel('Iceberg melt rate (m yr^{-1})','fontsize',20);
            set(gca,'xlim',[0 800],'xtick',[0:200:800],'xticklabel',[0:200:800]);
            set(ylb,'position',[-100 350 -1]);
        end
        
        %adjust subplot positions
        if j == 1
            if mod(plot_loc(i),2) ~= 0
                set(gca,'position',[plotpos(1)-0.03 plotpos(2)+0.04/(plot_loc(i)+1) 1.1*plotpos(3) 1.1*plotpos(4)]);
            else
                set(gca,'position',[plotpos(1) plotpos(2)+0.04/(plot_loc(i)) 1.1*plotpos(3) 1.1*plotpos(4)]);
            end
        end
        clear plotpos;
        if plot_loc(i) == 14
            if j == length(meltinfo)
                legpos = get(legsub,'position'); set(legsub,'position',[0.56 0.11 legpos(3) legpos(4)]); clear legpos;
            end
        end
        
        %individual region plots
        figure(figureE);
        %color-coded according to region w/ edges color-coded by date
%         errorbar(draft,dVdt./Asub,abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,[marker(i),'k'],'markerfacecolor',plot_cmap(i,:),'markeredgecolor',plot_marker(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markersize',symbol_size); hold on;
%         indpl(j) = plot(draft,dVdt./Asub,[marker(i),'k'],'markerfacecolor',plot_cmap(i,:),'markeredgecolor',plot_marker(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markersize',symbol_size,'linewidth',1.5); hold on;
        %plot face color color-coded by date
        errorbar(draft,365*dVdt./Asub,365*abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),365*abs(dVdt./Asub).*sqrt((dVdt_uncert./dVdt).^2 + (Asub_uncert./Asub).^2),draft_uncert,draft_uncert,[marker(i),'k'],'markerfacecolor',plot_marker(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markeredgecolor','k','markersize',symbol_size); hold on;
        indpl(j) = plot(draft,365*dVdt./Asub,[marker(i),'k'],'markerfacecolor',plot_marker(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markeredgecolor','k','markersize',symbol_size); hold on;
%         indpl(j) = plot(draft,dVdt./Asub,[marker(i),'k'],'markerfacecolor',plot_marker(round(nanmean(plot_yrs))-(min(yrs)-1),:),'markeredgecolor','k','markersize',symbol_size,'linewidth',1.5); hold on;
        daterange(j,:) = [meltinfo(j).name(end-49:end-46),'/',meltinfo(j).name(end-45:end-44),'/',meltinfo(j).name(end-43:end-42),'-',meltinfo(j).name(end-34:end-31),'/',meltinfo(j).name(end-30:end-29),'/',meltinfo(j).name(end-28:end-27)];
        if j==length(meltinfo)
            xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
            if max(xlims) > 500
                set(gca,'xlim',[0 max(xlims)],'xtick',[0:100:round(max(xlims))],'xticklabel',[0:100:round(max(xlims))],...
                    'ylim',[0 max(ylims)],'ytick',[0:round(max(ylims)/3):max(ylims)],'yticklabel',[0:round(max(ylims)/3):max(ylims)],'fontsize',20); grid on;
            else
                set(gca,'xlim',[0 max(xlims)],'xtick',[0:50:round(max(xlims))],'xticklabel',[0:50:round(max(xlims))],...
                    'ylim',[0 max(ylims)],'ytick',[0:round(max(ylims)/3):max(ylims)],'yticklabel',[0:round(max(ylims)/3):max(ylims)],'fontsize',20); grid on;
            end
            xlabel('Draft (m b.s.l.)','fontsize',20); ylabel('Melt rate (m yr^{-1})','fontsize',20);
            leg = legend(indpl,daterange); set(leg,'fontsize',16,'location','northwest');
            saveas(figureE,[char(region(i)),'_iceberg-meltrate-v-draft.eps'],'epsc'); saveas(figureE,[char(region(i)),'_iceberg-meltrate-v-draft.png'],'png');
        end
        clear dt xo yo zo rhoo Vo xf yf zf rhof Vf coregzo coregzf dz dz_sigma dVdt dVdt_uncert draft draft_uncert Asurf Asurf_uncert Asub Asub_uncert;
    end
    figure(figureC);
    text(0.6*max(get(gca,'xlim')),0.85*max(get(gca,'ylim')),char(disp_names(i)),'fontsize',16);
    %add a map-view figure
    if ~isempty(landsats)
        [A,S] = geotiffread(landsats(1).name);
        im.x = S.XWorldLimits(1):S.SampleSpacingInWorldX:S.XWorldLimits(2);
        im.y = S.YWorldLimits(2):-S.SampleSpacingInWorldY:S.YWorldLimits(1);
        im.z = double(A); clear A S;
        figureD = figure; set(gcf,'position',[650 150 900 500]);
        imagesc(im.x,im.y,im.z); axis xy equal; colormap(gray(10001)); hold on;
        for j = 1:length(meltrate)
            symbol_color = ceil(2*(meltrate(j))*365);
%             if i>=8
                if symbol_color > length(highmelt_cmap); symbol_color = length(highmelt_cmap); end
                plot(nanmean([xcoord_o(j) xcoord_f(j)]),nanmean([ycoord_o(j) ycoord_f(j)]),[marker(i),'k'],'markerfacecolor',highmelt_cmap(symbol_color,:),'markersize',round(keeld(j)/30)+10); hold on;
%             else
%                 if symbol_color > length(lowmelt_cmap); symbol_color = length(lowmelt_cmap); end
%                 plot(nanmean([xcoord_o(j) xcoord_f(j)]),nanmean([ycoord_o(j) ycoord_f(j)]),'ok','markerfacecolor',lowmelt_cmap(symbol_color,:),'markersize',round(keeld(j)/30)+10); hold on;
%             end
            clear symbol_color;
        end
        if max(xcoord_o)-min(xcoord_o)+10000 < 50000
            set(gca,'xlim',[min(xcoord_o)-5000 max(xcoord_o)+5000],'xtick',[(ceil(min(xcoord_o)/1000)*1000-5000):5000:(floor(max(xcoord_o)/1000)*1000+5000)],'xticklabel',[(ceil(min(xcoord_o)/1000)-5):5:(floor(max(xcoord_o)/1000)+5)],...
                'ylim',[min(ycoord_o)-5000 max(ycoord_o)+5000],'ytick',[(ceil(min(ycoord_o)/1000)*1000-5000):5000:(floor(max(ycoord_o)/1000)*1000+5000)],'yticklabel',[(ceil(min(ycoord_o)/1000)-5):5:(floor(max(ycoord_o)/1000)+5)],...
                'fontsize',20);
        else
            set(gca,'xlim',[min(xcoord_o)-5000 max(xcoord_o)+5000],'xtick',[(ceil(min(xcoord_o)/1000)*1000-5000):10000:(floor(max(xcoord_o)/1000)*1000+5000)],'xticklabel',[(ceil(min(xcoord_o)/1000)-5):10:(floor(max(xcoord_o)/1000)+5)],...
                'ylim',[min(ycoord_o)-5000 max(ycoord_o)+5000],'ytick',[(ceil(min(ycoord_o)/1000)*1000-5000):5000:(floor(max(ycoord_o)/1000)*1000+5000)],'yticklabel',[(ceil(min(ycoord_o)/1000)-5):5:(floor(max(ycoord_o)/1000)+5)],...
                'fontsize',20);
        end
        xlabel('Easting (km)','fontsize',20); ylabel('Northing (km)','fontsize',20);
        %resize vertical dimension to maximize figure window usage
        xlims = get(gca,'xlim'); ylims = get(gca,'ylim'); figpos = get(gcf,'position');
        set(gcf,'position',[figpos(1) figpos(2) figpos(3) (max(ylims)-min(ylims))/(max(xlims)-min(xlims))*figpos(3)]);
        %add legends
        if ~strmatch('Thwaites',char(region(i)))
            rectangle('position',[min(xlims)+0.25*(max(xlims)-min(xlims)) min(ylims)+0.025*(max(ylims)-min(ylims)) 0.15*(max(xlims)-min(xlims)) 0.175*(max(ylims)-min(ylims))],'curvature',[0,0],'facecolor','w','linewidth',2);
            plot(min(xlims)+0.275*(max(xlims)-min(xlims)),min(ylims)+0.17*(max(ylims)-min(ylims)),[marker(i),'k'],'markersize',round(50/30+10),'markerfacecolor','w'); text(min(xlims)+0.30*(max(xlims)-min(xlims)),min(ylims)+0.1675*(max(ylims)-min(ylims)),'50 m','fontsize',16);
            plot(min(xlims)+0.275*(max(xlims)-min(xlims)),min(ylims)+0.12*(max(ylims)-min(ylims)),[marker(i),'k'],'markersize',round(150/30+10),'markerfacecolor','w'); text(min(xlims)+0.30*(max(xlims)-min(xlims)),min(ylims)+0.1175*(max(ylims)-min(ylims)),'150 m','fontsize',16);
            plot(min(xlims)+0.275*(max(xlims)-min(xlims)),min(ylims)+0.06*(max(ylims)-min(ylims)),[marker(i),'k'],'markersize',round(300/30+10),'markerfacecolor','w'); text(min(xlims)+0.30*(max(xlims)-min(xlims)),min(ylims)+0.0575*(max(ylims)-min(ylims)),'300 m','fontsize',16);
            rectangle('position',[min(xlims)+0.025*(max(xlims)-min(xlims)) min(ylims)+0.025*(max(ylims)-min(ylims)) 0.20*(max(xlims)-min(xlims)) 0.25*(max(ylims)-min(ylims))],'curvature',[0,0],'facecolor','w','linewidth',2);
            %             symbol_color = round(([0.001 0.01 0.1])*1000)+1; symbol_color(symbol_color>length(highmelt_cmap)) = length(highmelt_cmap);
            for k = 1:length(highmelt_cmap)
                plot([min(xlims)+0.05*(max(xlims)-min(xlims)) min(xlims)+0.09*(max(xlims)-min(xlims))],...
                    [min(ylims)+0.245*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap)) min(ylims)+0.245*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap))],'-','linewidth',1.5,'color',highmelt_cmap(k,:));
            end
            text(min(xlims)+0.10*(max(xlims)-min(xlims)),min(ylims)+0.25*(max(ylims)-min(ylims))-((0.004*(max(ylims)-min(ylims)))),'1 m yr^{-1}','fontsize',16);
            text(min(xlims)+0.10*(max(xlims)-min(xlims)),min(ylims)+0.25*(max(ylims)-min(ylims))-((0.04*(max(ylims)-min(ylims)))),'10 m yr^{-1}','fontsize',16);
            text(min(xlims)+0.10*(max(xlims)-min(xlims)),min(ylims)+0.25*(max(ylims)-min(ylims))-((0.20*(max(ylims)-min(ylims)))),'50 m yr^{-1}','fontsize',16);
        else
            rectangle('position',[min(xlims)+0.80*(max(xlims)-min(xlims)) min(ylims)+0.725*(max(ylims)-min(ylims)) 0.15*(max(xlims)-min(xlims)) 0.175*(max(ylims)-min(ylims))],'curvature',[0,0],'facecolor','w','linewidth',2);
            plot(min(xlims)+0.825*(max(xlims)-min(xlims)),min(ylims)+0.87*(max(ylims)-min(ylims)),[marker(i),'k'],'markersize',round(50/30+10),'markerfacecolor','w'); text(min(xlims)+0.85*(max(xlims)-min(xlims)),min(ylims)+0.8675*(max(ylims)-min(ylims)),'50 m','fontsize',16);
            plot(min(xlims)+0.825*(max(xlims)-min(xlims)),min(ylims)+0.82*(max(ylims)-min(ylims)),[marker(i),'k'],'markersize',round(150/30+10),'markerfacecolor','w'); text(min(xlims)+0.85*(max(xlims)-min(xlims)),min(ylims)+0.8175*(max(ylims)-min(ylims)),'150 m','fontsize',16);
            plot(min(xlims)+0.825*(max(xlims)-min(xlims)),min(ylims)+0.76*(max(ylims)-min(ylims)),[marker(i),'k'],'markersize',round(300/30+10),'markerfacecolor','w'); text(min(xlims)+0.85*(max(xlims)-min(xlims)),min(ylims)+0.7575*(max(ylims)-min(ylims)),'300 m','fontsize',16);
            rectangle('position',[min(xlims)+0.575*(max(xlims)-min(xlims)) min(ylims)+0.725*(max(ylims)-min(ylims)) 0.20*(max(xlims)-min(xlims)) 0.25*(max(ylims)-min(ylims))],'curvature',[0,0],'facecolor','w','linewidth',2);
            %             symbol_color = round(([0.001 0.01 0.1])*1000)+1; symbol_color(symbol_color>length(highmelt_cmap)) = length(highmelt_cmap);
            for k = 1:length(highmelt_cmap)
                plot([min(xlims)+0.60*(max(xlims)-min(xlims)) min(xlims)+0.64*(max(xlims)-min(xlims))],...
                    [min(ylims)+0.945*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap)) min(ylims)+0.945*(max(ylims)-min(ylims))-k*((0.20*(max(ylims)-min(ylims)))/length(highmelt_cmap))],'-','linewidth',1.5,'color',highmelt_cmap(k,:));
            end
            text(min(xlims)+0.65*(max(xlims)-min(xlims)),min(ylims)+0.95*(max(ylims)-min(ylims))-((0.004*(max(ylims)-min(ylims)))),'1 m yr^{-1}','fontsize',16);
            text(min(xlims)+0.65*(max(xlims)-min(xlims)),min(ylims)+0.95*(max(ylims)-min(ylims))-((0.04*(max(ylims)-min(ylims)))),'10 m yr^{-1}','fontsize',16);
            text(min(xlims)+0.65*(max(xlims)-min(xlims)),min(ylims)+0.95*(max(ylims)-min(ylims))-((0.20*(max(ylims)-min(ylims)))),'50 m yr^{-1}','fontsize',16);
        end
        clear im; clear symbol_color;
        saveas(figureD,[char(region(i)),'_melt-map.eps'],'epsc'); saveas(figureD,[char(region(i)),'_melt-map.png'],'png');
    end
    clear xlims ylims indpl daterange;
    
    %save all dates and coords to a tab-delimited text file
    cd ..
    date_coords = [date_o xcoord_o ycoord_o date_f xcoord_f ycoord_f];
    if i ==1
        dlmwrite('Antarctic-iceberg-PScoords.txt',date_coords,'delimiter','\t','precision',8);
    else
        dlmwrite('Antarctic-iceberg-PScoords.txt',date_coords,'-append','delimiter','\t','precision',8);
    end
    clear date_coords;
    clear subpl;
end
figure(figureA);
set(gca,'xlim',[-28e5 28e5],'xtick',[-32e5:8e5:32e5],'xticklabel',[-3200:800:3200],...
    'ylim',[-24e5 24e5],'ytick',[-24e5:8e5:24e5],'yticklabel',[-2400:800:2400],'fontsize',20); grid off;
xlabel('Easting (km)','fontsize',20); ylabel('Northing (km)','fontsize',20);
graticuleps(-50:-5:-90,-180:30:180);
text(0,6.5e5,['85',char(176),'S'],'fontsize',16); text(0,12.0e5,['80',char(176),'S'],'fontsize',16); 
text(0,17.5e5,['75',char(176),'S'],'fontsize',16); text(0,23.0e5,['70',char(176),'S'],'fontsize',16);
text(-16.5e5,25.25e5,['-30',char(176),'E'],'fontsize',16); text(12.5e5,25.25e5,['30',char(176),'E'],'fontsize',16); 
colormap(gca,im_cmap);
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
saveas(gcf,'Antarctic-iceberg-map.eps','epsc'); saveas(gcf,'Antarctic-iceberg-map.png','png');

%save the subplots containing all data
% Wfit = fit(WAsub,W365*dVdt,'poly1'); Efit = fit(EAsub,E365*dVdt,'poly1'); 
figure(figureB);
subplot(sub1b);
% plot([0 max([WAsub])],Wfit.p1.*[0 max([WAsub])] + Wfit.p2,'-k','linewidth',2); hold on;
% plot([0 max([EAsub])],Efit.p1.*[0 max([EAsub])] + Wfit.p2,'--k','linewidth',2); hold on;
leg1 = legend(pl,[char(leg_names)]); set(leg1,'location','westoutside','fontsize',16); set(leg1,'position',[0.03    0.2756    0.0685    0.4838]);
set(gca,'xlim',[0 7e6],'xtick',[0:1e6:7e6],'xticklabel',[0:1:7],...
    'ylim',[0 6.8],'ytick',[0:1:7],'yticklabel',[0:1:7],'fontsize',20); grid on;
xlabel('Submerged area (km^2)','fontsize',20); ylabel('Meltwater flux (m^3 s^{-1})','fontsize',20);
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
text(0.05*max(xlims),0.95*max(ylims),'a) ','color','k','fontsize',20);
set(sub1b,'position',[0.1600    0.1400    0.2134    0.8050]);
subplot(sub2b);
set(gca,'xlim',[0 800],'xtick',[0:200:800],'xticklabel',[0:200:800],...
    'ylim',[0 72],'ytick',[0:24:72],'yticklabel',[0:24:72],'fontsize',20); grid on;
xlabel('Draft (m b.s.l.)','fontsize',20); ylabel('Melt rate (m yr^{-1})','fontsize',20);
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
text(0.05*max(xlims),0.95*max(ylims),'b) ','color','k','fontsize',20);
rectangle('position',[0 0 600 24],'linewidth',2,'edgecolor','k','linestyle','--'); hold on;
set(sub2b,'position',[0.4408    0.1400    0.2134    0.8050]);
subplot(sub3b);
set(gca,'xlim',[0 800],'xtick',[0:100:800],'xticklabel',[0:100:800],...
    'ylim',[0 24],'ytick',[0:8:24],'yticklabel',[0:8:24],'fontsize',20); grid on;
xlabel('Draft (m b.s.l.)','fontsize',20); ylabel('Melt rate (m yr^{-1})','fontsize',20);
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
text(0.05*max(xlims),0.95*max(ylims),'c) ','color','k','fontsize',20);
% rectangle('position',[0 0 600 0.085],'linewidth',2,'edgecolor','k','linestyle','--'); hold on;
set(sub3b,'position',[0.7216    0.1400    0.2134    0.8050]);
saveas(gcf,'Antarctic-iceberg-lumped-plots.eps','epsc'); saveas(gcf,'Antarctic-iceberg-lumped-plots.png','png');

%save the subplots sorted by study site
figure(figureC);
% % leg2 = legend(yrpl,num2str(yrs')); set(leg2,'position',[0.92    0.35    0.045    0.26]);
saveas(gcf,'Antarctic-iceberg-subplots.eps','epsc'); saveas(gcf,'Antarctic-iceberg-subplots.png','png');


% %% replot Larsen A & B to hone-in on potential variations btw dates
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
%                 'ylim',[0 0.04],'ytick',[0:0.01:0.04],'yticklabel',[0:1:4],'fontsize',20); grid on;
%     text(0.1*max(get(gca,'xlim')),0.9*max(get(gca,'ylim')),char(region(i)),'fontsize',20);
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
% %             'fontsize',20);
% %         xlabel('Easting (km)','fontsize',20); ylabel('Northing (km)','fontsize',20);
% %         clear im;
% %         saveas(figureD,[char(region(i)),'_melt-map.eps'],'epsc'); saveas(figureD,[char(region(i)),'_melt-map.png'],'png');
% %     end
%     cd ..
% end