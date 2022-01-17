%%% Compare Iceberg melt rates to Parameterized melt rates

%%inputs:
%1) Antarctic-icebergmelt-comparison.mat: melt structure w/ iceberg data & nearby CTD and/or APB data
%2) SOSE data: bsose_i134_2013to2019_monthly_to60S_Salt.nc (salinity; psu), bsose_i134_2013to2019_monthly_to60S_Theta.nc (temp; degrees C),
%bsose_i134_2013to2019_monthly_to60S_Uvel.nc (x velocity; m/s), & bsose_i134_2013to2019_monthly_to60S_Vvel.nc (y velocity; m/s)

%%outputs:


%%dependencies:
%1) ps2wgs.m
%2) cmocean.m

%%%
%% initiate
clearvars; close all;
addpath('/users/ellynenderlin/mfiles','/users/ellynenderlin/mfiles/general','/users/ellynenderlin/mfiles/general/cmocean');

%set variables
rho_sw = 1026; %density of water (kg/m^3) 
diffM = 1.00*10^-2/(100^2); %diffusivity of momentum (m^2/s)
diffT = 1.42*10^-3/(100^2); %diffusivity of temperature (m^2/s)
c = 4.182*1000; %heat capacity (J/(kg*C))
L = 334*1000; %latent heat of melting (J/kg)

%specify plot params
% marker = ['d','d','d','d','d','d','d','s','s','s','s','s','s','s','s']; %set-up the different marker styles for E and W
marker = ['s','s','s','s','s','s','s','s','s','s','s','s','s','s','s']; %set-up the different marker styles for E and W
plot_names = [{'i)'},{'j)'},{'k)'},{'l)'},{'m)'},{'n)'},{'o)'},{'h)'},{'g)'},{'f)'},{'e)'},{'d)'},{'c)'},{'b)'},{'a)'}]; %plot letters
plot_loc = [2,4,6,8,10,12,14,15,13,11,9,7,5,3,1];

%load the iceberg & ocean observation data
cd /Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt
load Antarctic-icebergmelt-comparison.mat; 


%% read-in SOSE data (T,S,U&V) & observed ocean T&S

%load the SOSE ocean model outputs
disp('Load SOSE data... (this takes a while)')
cd /Users/ellynenderlin/Research/miscellaneous/SOSE/
SOSEt = ncread('bsose_i134_2013to2019_monthly_to60S_Salt.nc','time');
SOSEx = ncread('bsose_i134_2013to2019_monthly_to60S_Salt.nc','XC');
SOSEy = ncread('bsose_i134_2013to2019_monthly_to60S_Salt.nc','YC'); SOSEy = SOSEy';
SOSEz = ncread('bsose_i134_2013to2019_monthly_to60S_Salt.nc','Z'); SOSEd = -SOSEz; clear SOSEz;
SOSEbottom = ncread('bsose_i134_2013to2019_monthly_to60S_Salt.nc','Depth'); %maximum depth (m)
SOSES = ncread('bsose_i134_2013to2019_monthly_to60S_Salt.nc','SALT');
SOSET = ncread('bsose_i134_2013to2019_monthly_to60S_Theta.nc','THETA');
SOSEU = ncread('bsose_i134_2013to2019_monthly_to60S_Uvel.nc','UVEL');
SOSEV = ncread('bsose_i134_2013to2019_monthly_to60S_Vvel.nc','VVEL');
% SOSEdate = 2013+single(SOSEt)./31536000; clear SOSEt;
SOSEdate = [2013+[0:1:11]/12 2014+[0:1:11]/12 2015+[0:1:11]/12 2016+[0:1:11]/12 2017+[0:1:11]/12 2018+[0:1:11]/12 2019+[0:1:11]/12];
SOSExgrid = repmat(SOSEx,1,length(SOSEy)); SOSEygrid = repmat(SOSEy,length(SOSEx),1); 
disp('Read-in SOSE data');

%create time-averaged SOSE datasets to use when iceberg observations fall
%outside the 2013-2019 SOSE date range
for k = 1:12
    SOSESmean(:,:,:,k) = squeeze(nanmean(SOSES(:,:,:,k:12:end),4)); SOSETmean(:,:,:,k) = squeeze(nanmean(SOSET(:,:,:,k:12:end),4));
    SOSEUmean(:,:,:,k) = squeeze(nanmean(SOSEU(:,:,:,k:12:end),4)); SOSEVmean(:,:,:,k) = squeeze(nanmean(SOSEV(:,:,:,k:12:end),4));
end
SOSEtmean = [0:1:11]/12;
disp('Calculated SOSE monthly means from 2013 through 2019');

%create an avg temp map 
Tmap = squeeze(nanmean(nanmean(SOSET,4),3)); Smap = squeeze(nanmean(nanmean(SOSES,4),3));
Velmap = squeeze(nanmean(nanmean(SOSEU,4),3));
Tmap(Velmap==0) = NaN;

%assume the iceberg is moving at the same speed as the water averaged
%over some depth from the surface
depthref = find(SOSEd<=50,1,'last');

%calculate the thickness & depth range of each SOSE vertical layer
SOSEH(1) = 2*SOSEd(1); SOSEdbottom(1) = 2*SOSEd(1); 
for j=2:length(SOSEd)
    SOSEH(j) = 2*(SOSEd(j)-SOSEdbottom(j-1));
    SOSEdbottom(j) = SOSEdbottom(j-1)+SOSEH(j);
end
SOSEH = SOSEH';

%% estimate melt using SOSE data
close all; drawnow;
cd /Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt

%set-up figure for subplots
fig1 = figure; set(gcf,'position',[50 50 800 1000]); 
yr_cmap = cmocean('dense',9); yr_cmap = yr_cmap(2:end,:);
Temp_cmap = cmocean('thermal',600); Temp_cmap(1,:) = [1 1 1];
imagesc(Tmap); colormap(gca,Temp_cmap); hold on;
fig2 = figure; set(gcf,'position',[350 50 800 1000]);
fig3 = figure; set(gcf,'position',[650 50 800 1000]);

%create a textfile w/ headers to export data
combined_filename = fullfile(pwd,'Antarctic-icebergs-summary-data.txt');
fid = fopen(combined_filename,'wt');
fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','x (PS E)','y (PS N)','first year','last year','median draft (m)','MAD draft (m)',...
    'iceberg temp. (C)','median density (kg/m^3)','MAD density (kg/m^3)','median obs. meltrate (m/yr)','MAD obs. meltrate (m/yr)','median param. meltrate (m/yr)','MAD param. meltrate (m/yr)',...
    'median HJ-WC param. meltrates (m/yr)','median obs.-SOSE temp. (C)','median Urelative (m/s)','MAD Urelative (m/s)');
fclose(fid);

%estimate melt rates for each iceberg using SOSE data
disp('Calculating melt using SOSE data plugged into two melt parameterizations...');
for i = 1:length(melt)
    disp(char(melt(i).dispname)); 
    site_fig = figure; set(gcf,'position',[850 50 1600 400]);
    figure(fig3); subpl = subplot(8,2,plot_loc(i));

    %set-up an empty matrix to fill with median +/- MAD melt rates for each date
    meltrates = [];
    obs_dates = [melt(i).to melt(i).tf];
    [unique_dates,unique_refs,unique_valrefs] = unique(obs_dates,'rows');
    
    %loop through each iceberg
    for j = 1:length(melt(i).d)
        
        %identify the nearest SOSE cell with data
        [lon,lat] = ps2wgs(melt(i).x(j),melt(i).y(j),'StandardParallel',-71,'StandardMeridian',0); 
        lon(lon<0) = lon(lon<0)+360;
%         figure(fig1); plot(lat,lon,'xr','linewidth',2,'markersize',20); hold on;
        diffmap = sqrt((lon-SOSExgrid).^2 + (lat-SOSEygrid).^2);
        diffmap(Velmap==0 | isnan(Velmap)) = NaN;
        map_ind = find(diffmap == min(min(diffmap)));
        [x_ind,y_ind] = ind2sub(size(Tmap),map_ind);
        figure(fig1); plot(y_ind,x_ind,'+m','linewidth',2,'markersize',20); hold on;
        clear lon lat diffmap map_ind;
        
        %find the SOSE indices for the appropriate date range
        if melt(i).to(j) > 2013 && melt(i).tf(j) < 2020
            earlyref = find(abs(melt(i).to(j)-SOSEdate) == min(abs(melt(i).to(j)-SOSEdate)));
            lateref = find(abs(melt(i).tf(j)-SOSEdate) == min(abs(melt(i).tf(j)-SOSEdate)));
            Tprofs = squeeze(SOSET(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:lateref));
            Sprofs = squeeze(SOSES(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:lateref));
            Uprofs = squeeze(SOSEU(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:lateref));
            Vprofs = squeeze(SOSEV(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:lateref));
            
            clear earlyref lateref;
            
        elseif melt(i).to(j) < 2013
            %identify the partial years before 2013
            partyrs = ceil(12*(ceil(melt(i).to(j)) - melt(i).to(j)));
            
            %determine if the entire date range is before 2013
            if melt(i).tf(j) < 2013
                %identify the number of full years before 2013
                fullyrs = floor(melt(i).tf(j))-ceil(melt(i).to(j));
            
                %identify the end reference
                lateref = ceil(12*(melt(i).tf(j) - floor(melt(i).tf(j))));

                %combine averaged data for the full time period
                Tprofs = horzcat(squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                Sprofs = horzcat(squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                Uprofs = horzcat(squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                Vprofs = horzcat(squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                
            else
                %identify the number of full years before 2013
                fullyrs = 2013-ceil(melt(i).to(j));
                
                %identify the end reference
                lateref = find(abs(melt(i).tf(j)-SOSEdate) == min(abs(melt(i).tf(j)-SOSEdate)));
                
                %combine averaged data from before 2013 & post-2013 data
                Tprofs = horzcat(squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSET(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                Sprofs = horzcat(squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSES(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                Uprofs = horzcat(squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEU(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                Vprofs = horzcat(squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEV(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                
            end
            clear partyrs fullyrs lateref;
            
        elseif melt(i).tf(j) > 2019
            %identify the partial years after 2019
            partyrs = ceil(12*(melt(i).tf(j) - floor(melt(i).tf(j))));
            
            %determine if the entire date range is after 2019
            if melt(i).to(j) > 2019
                %identify the number of full years after 2019
                fullyrs = floor(melt(i).tf(j))-ceil(melt(i).to(j));
            
                %identify the start reference
                earlyref = ceil(12*(ceil(melt(i).to(j)) - melt(i).to(j)));

                %combine averaged data for the full time period
                Tprofs = horzcat(squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-earlyref:12)),...
                    fullyrs*squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                Sprofs = horzcat(squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-earlyref:12)),...
                    fullyrs*squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                Uprofs = horzcat(squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-earlyref:12)),...
                    fullyrs*squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                Vprofs = horzcat(squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-earlyref:12)),...
                    fullyrs*squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                
            else
                %identify the number of full years after 2019
                fullyrs = floor(melt(i).tf(j))-2019;
                
                %identify the start reference
                earlyref = find(abs(melt(i).to(j)-SOSEdate) == min(abs(melt(i).to(j)-SOSEdate)));
                
                %combine averaged data from before 2013 & post-2013 data
                Tprofs = horzcat(squeeze(SOSET(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:end)),...
                    fullyrs*squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                Sprofs = horzcat(squeeze(SOSES(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:end)),...
                    fullyrs*squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                Uprofs = horzcat(squeeze(SOSEU(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:end)),...
                    fullyrs*squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                Vprofs = horzcat(squeeze(SOSEV(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:end)),...
                    fullyrs*squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                
            end
            clear partyrs fullyrs earlyref;
            
        end
        
        %average the profiles over time
        Tavg=nanmean(Tprofs,2); Savg=nanmean(Sprofs,2); Uavg=nanmean(Uprofs,2); Vavg=nanmean(Vprofs,2); davg = SOSEd(1:length(Tavg));
        Tavg(Savg==0) = []; Uavg(Savg==0) = []; Vavg(Savg==0) = []; davg(Savg==0) = []; Savg(Savg==0) = []; 
        clear *profs;
        
        %determine the draft range that gives iceberg velocities that
        %minimize the difference between observed and modeled melt rates
        for k = 1:length(Tavg)
            depthref = k;
            
            %estimate iceberg velocity as the depth-averaged water velocity
            U_weight = nansum(Uavg(1:depthref).*SOSEH(1:depthref))./nansum(SOSEH(1:depthref));
            V_weight = nansum(Vavg(1:depthref).*SOSEH(1:depthref))./nansum(SOSEH(1:depthref));

            %estimate melt rate using the Holland & Jenkins/Silva melt equations
            [Mhj, T_sh, T_fp,a] = melt_forcedwater(Tavg,Savg,davg, sqrt((Uavg-U_weight).^2 + (Vavg-V_weight).^2), abs(melt(i).bergT));
            %calculate the average melt rate, weighting according to lateral & basal areas
            lat_area = melt(i).Asub(j)-melt(i).Asurf(j); bas_area = melt(i).Asurf(j);
            Mhj_weighted = nansum(Mhj.*SOSEH(1:length(Tavg)))./nansum(SOSEH(1:length(Tavg)));
            mHJ(k) = 31536000*((Mhj_weighted*lat_area)+(Mhj(end)*bas_area))./(lat_area+bas_area); %area-averaged melt rate (m yr^{-1})
            %estimate melt rate using the Weeks & Campbell melt equation
            Mwc = 0.037*((rho_sw/melt(i).rho(j))*diffM^(-7/15)*diffT^(2/3)*(c/L))*((sqrt((Uavg-U_weight).^2 + (Vavg-V_weight).^2).^(0.8).*(Tavg-T_fp))./((2*sqrt(melt(i).Asurf(j)/pi)).^0.2));
            Mwc_weighted = nansum(Mwc.*SOSEH(1:length(Tavg)))./nansum(SOSEH(1:length(Tavg)));
            mWC(k) = 31536000*((Mwc_weighted*lat_area)+(Mwc(end)*bas_area))./(lat_area+bas_area); %area-averaged melt rate (m yr^{-1})
            clear Mhj* Mwc* U_weight V_weight;
        end
        depthref = find(melt(i).m(j)./nanmean([mWC; mHJ]) == min(melt(i).m(j)./nanmean([mWC; mHJ])));
        melt(i).SOSE_icebergUref(j) = depthref; melt(i).SOSE_icebergUdepth(j) = SOSEdbottom(depthref);
        clear mWC mHJ;
        %estimate iceberg velocity as the depth-averaged water velocity
        U_weighted(j) = nansum(Uavg(1:depthref).*SOSEH(1:depthref))./nansum(SOSEH(1:depthref));
        V_weighted(j) = nansum(Vavg(1:depthref).*SOSEH(1:depthref))./nansum(SOSEH(1:depthref));
        %estimate melt rate using the Holland & Jenkins/Silva melt equations
        [Mhj, T_sh, T_fp,a] = melt_forcedwater(Tavg,Savg,davg, sqrt((Uavg-U_weighted(j)).^2 + (Vavg-V_weighted(j)).^2), abs(melt(i).bergT));
        %calculate the average melt rate, weighting according to lateral & basal areas
        lat_area = melt(i).Asub(j)-melt(i).Asurf(j); bas_area = melt(i).Asurf(j);
        Mhj_weighted = nansum(Mhj.*SOSEH(1:length(Tavg)))./nansum(SOSEH(1:length(Tavg)));
        melt(i).mHJ(j) = 31536000*((Mhj_weighted*lat_area)+(Mhj(end)*bas_area))./(lat_area+bas_area); %area-averaged melt rate (m yr^{-1})
        %estimate melt rate using the Weeks & Campbell melt equation
        Mwc = 0.037*((rho_sw/melt(i).rho(j))*diffM^(-7/15)*diffT^(2/3)*(c/L))*((sqrt((Uavg-U_weighted(j)).^2 + (Vavg-V_weighted(j)).^2).^(0.8).*(Tavg-T_fp))./((2*sqrt(melt(i).Asurf(j)/pi)).^0.2));
        Mwc_weighted = nansum(Mwc.*SOSEH(1:length(Tavg)))./nansum(SOSEH(1:length(Tavg)));
        melt(i).mWC(j) = 31536000*((Mwc_weighted*lat_area)+(Mwc(end)*bas_area))./(lat_area+bas_area); %area-averaged melt rate (m yr^{-1})
        clear Mhj* Mwc*;
        
        
        %melt rates using observed T&S but SOSE vels
        if ~isempty(melt(i).oceanTavg)
            %only have depth-averaged T&S saved: create T&S avg profiles
            Tobs = repmat(melt(i).oceanTavg(j),size(Uavg)); Sobs = repmat(melt(i).oceanSavg(j),size(Uavg)); 
            %create a glacier length profile
            Lobs = repmat(2*sqrt(melt(i).Asurf(j)/pi),size(Uavg));
            
            %estimate melt rate using the Holland & Jenkins/Silva melt equations
            [Mhj,~,~,~] = melt_forcedwater(Tobs,Sobs,davg, sqrt((Uavg-U_weighted(j)).^2 + (Vavg-V_weighted(j)).^2), abs(melt(i).bergT));
            %calculate the average melt rate, weighting according to lateral & basal areas
            Mhj_weighted = nansum(Mhj.*SOSEH(1:length(Tavg)))./nansum(SOSEH(1:length(Tavg)));
            melt(i).mHJ_TSobs(j) = 31536000*((Mhj_weighted*lat_area)+(Mhj(end)*bas_area))./(lat_area+bas_area); %area-averaged melt rate (m yr^{-1})
            %estimate melt rate using the Weeks & Campbell melt equation
%             T_fp = a(1).*melt(i).oceanSavg(j) + a(2) + a(3).*melt(i).d(j)/2; % fp temp
            Mwc = 0.037*((rho_sw/melt(i).rho(j))*diffM^(-7/15)*diffT^(2/3)*(c/L))*((sqrt((Uavg-U_weighted(j)).^2 + (Vavg-V_weighted(j)).^2).^(4/5).*(Tobs-T_fp))./((Lobs).^(1/5)));
            Mwc_weighted = nansum(Mwc.*SOSEH(1:length(Tavg)))./nansum(SOSEH(1:length(Tavg)));
            melt(i).mWC_TSobs(j) = 31536000*((Mwc_weighted*lat_area)+(Mwc(end)*bas_area))./(lat_area+bas_area); %area-averaged melt rate (m yr^{-1})
            
            %calculate the relative velocity needed to match observed melt rates
            melt(i).oceanUWC_est(j) = abs((((melt(i).m(j)/86400)*((Lobs(1)).^(1/5)))/(0.037*((rho_sw/melt(i).rho(j))*diffM^(-7/15)*diffT^(2/3)*(c/L))*(melt(i).oceanTavg(j)-nanmean(T_fp))))^(5/4)); 
            melt(i).oceanUHJ_est(j) = abs(melt(i).mHJ_TSobs(j)/melt(i).mWC_TSobs(j))*abs(melt(i).oceanUWC_est(j))^(4/5);
            clear *obs Mwc* Mhj*;
            
            %report the ratio of depth-averaged observed:SOSE temp
            temp_weighted(j) = nansum(Tavg.*SOSEH(1:length(Tavg)))./nansum(SOSEH(1:length(Tavg)));
        end
        clear T_fp;
        
        %plot the mps vs observed melt rate
        figure(fig3); subplot(subpl);
%         plT(1) = plot(nanmean(Tavg),365*melt(i).m(j),'s','markerfacecolor','none','color',nanmean([33,113,181; 107,174,214])/255,'linewidth',1.5); hold on;
%         plT(2) = plot(Tavg(end),365*melt(i).m(j),'^','markerfacecolor','none','color',nanmean([33,113,181; 107,174,214])/255,'linewidth',1.5); hold on;
        plT(1) = plot(((nanmean(Tavg)*lat_area)+(Tavg(end)*bas_area))./(lat_area+bas_area),365*melt(i).m(j),'^','markerfacecolor','none','color',[0.25 0.25 0.25],'linewidth',1.5); hold on;
        
        %create date-specific subplots for each glacier
        figure(site_fig); 
        subplot(1,3,1); set(gca,'fontsize',20); grid on; %observed melt rate vs draft
        if j == 1
            for k = 1:8
                plyr(k) = plot(melt(i).d(j),365*melt(i).m(j),[marker(i),'k'],'markerfacecolor',yr_cmap(k,:),'markersize',12); hold on;
            end
        end
%         errorbar(melt(i).d(j),365*melt(i).m(j),melt(i).d_uncert(j),melt(i).d_uncert(j),...
%             365*melt(i).m_uncert(j),365*melt(i).m_uncert(j),'sk'); hold on;
        plot(melt(i).d(j),365*melt(i).m(j),[marker(i),'k'],'markerfacecolor',yr_cmap(round(nanmean([melt(i).to(j) melt(i).tf(j)]))-2012,:),'markersize',12); hold on;
        set(gca,'ylim',[0 ceil(max(365*melt(i).m))],'xlim',[0 ceil(max(melt(i).d)+0.05*max(melt(i).d))]);
        ylabel('Iceberg melt rate (m yr^{-1})','fontsize',20);
        xlabel('Iceberg draft (m b.s.l.)','fontsize',20);
        subplot(1,3,2); set(gca,'fontsize',20); grid on; %Parameterized melt rate using all SOSE data & SOSE w/ observed temps 
        plot(melt(i).d(j),nanmean([melt(i).mWC(j); melt(i).mHJ(j)]),'^k','markerfacecolor',yr_cmap(round(nanmean([melt(i).to(j) melt(i).tf(j)]))-2012,:),'markersize',12); hold on;
        if ~isempty(melt(i).oceanTavg)
            plot(melt(i).d(j),nanmean([melt(i).mWC_TSobs(j); melt(i).mHJ_TSobs(j)]),'x','color',yr_cmap(round(nanmean([melt(i).to(j) melt(i).tf(j)]))-2012,:),'markersize',14,'linewidth',1.5); hold on;
        end
        set(gca,'ylim',[0 ceil(max(365*melt(i).m))],'xlim',[0 ceil(max(melt(i).d)+0.05*max(melt(i).d))]);
        ylabel('Parameterized iceberg melt rate (m yr^{-1})','fontsize',20);
        xlabel('Iceberg draft (m b.s.l.)','fontsize',20);
        subplot(1,3,3); set(gca,'fontsize',20); grid on;
        vel_weighted(j) = nansum(sqrt((Uavg-U_weighted(j)).^2 + (Vavg-V_weighted(j)).^2).*SOSEH(1:length(Tavg)))./nansum(SOSEH(1:length(Tavg)));
        plot(melt(i).d(j),vel_weighted(j),'^k','markerfacecolor',yr_cmap(round(nanmean([melt(i).to(j) melt(i).tf(j)]))-2012,:),'markersize',12); hold on;
        if ~isempty(melt(i).oceanTavg)
%             errorbar(melt(i).d(j),365*melt(i).m(j),...
%                 nanmean([melt(i).oceanUHJ_est(j) melt(i).oceanUWC_est(j)])-min([melt(i).oceanUHJ_est(j) melt(i).oceanUWC_est(j)]),max([melt(i).oceanUHJ_est(j) melt(i).oceanUWC_est(j)])-nanmean([melt(i).oceanUHJ_est(j) melt(i).oceanUWC_est(j)]),...
%                 365*melt(i).m_uncert(j),365*melt(i).m_uncert(j),'sk'); hold on;
            plot(melt(i).d(j),nanmean([melt(i).oceanUHJ_est(j) melt(i).oceanUWC_est(j)]),'x','color',yr_cmap(round(nanmean([melt(i).to(j) melt(i).tf(j)]))-2012,:),'markersize',14,'linewidth',1.5); hold on;
        end
        set(gca,'xlim',[0 ceil(max(melt(i).d)+0.05*max(melt(i).d))]);
        xlabel('Iceberg draft (m b.s.l.)','fontsize',20);
        ylabel('Relative velocity (m s^{-1})','fontsize',20);
        clear lat_area bas_area;
        
        
        clear *_ind *profs *avg M* T_*;
    end
    drawnow;
    
    %save meltrates to structure
    for k = 1:length(unique_refs)
        %dates
        meltdates(k,:) = unique_dates(k,:);
        %estimated draft
        meltrates(k,1) = nanmedian(melt(i).d(unique_valrefs==k));
        meltrates(k,2) = mad(melt(i).d(unique_valrefs==k),1);
        %observed melt rates
        meltrates(k,3) = 365*nanmedian(melt(i).m(unique_valrefs==k));
        meltrates(k,4) = 365*mad(melt(i).m(unique_valrefs==k),1);
        %parameterized melt rates with SOSE
        meltrates(k,5) = nanmedian(nanmean([melt(i).mHJ(unique_valrefs==k); melt(i).mWC(unique_valrefs==k)],1));
        meltrates(k,6) = mad(nanmean([melt(i).mHJ(unique_valrefs==k); melt(i).mWC(unique_valrefs==k)],1),1);
        
    end
    %export date-specific, site-specific information as a text file
    filename = fullfile(pwd,[char(melt(i).dispname),'-observed-parameterized-melt.txt']);
    fid = fopen(filename,'wt');
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n','time1','time2','median draft (m)',...
        'MAD draft (m)','median obs. meltrate (m/yr)','MAD obs. meltrate (m/yr)','median param. meltrate (m/yr)','MAD param. meltrate (m/yr)');
    fclose(fid);
    dlmwrite(filename,horzcat(meltdates, meltrates),'delimiter','\t',...
        'precision','%.2f','-append');
    clear meltdates meltrates;
    
    %compile & export time-averaged data to the full dataset textfile
    meltdata(1,1) = nanmedian(melt(i).x(~isnan(melt(i).m))); %median easting
    meltdata(1,2) = nanmedian(melt(i).y(~isnan(melt(i).m))); %median northing
    meltdata(1,3:4) = [min(melt(i).to) max(melt(i).tf)]; %date range
    meltdata(1,5) = nanmedian(melt(i).d); %median draft
    meltdata(1,6) = mad(melt(i).d(~isnan(melt(i).m)),1); %mad draft
    meltdata(1,7) = melt(i).bergT; %iceberg temp
    meltdata(1,8) = nanmedian(melt(i).rho); %median iceberg density
    meltdata(1,9) = mad(melt(i).rho(~isnan(melt(i).m)),1); %mad iceberg density
    meltdata(1,10) = 365*nanmedian(melt(i).m); %median observed melt rate
    meltdata(1,11) = 365*mad(melt(i).m(~isnan(melt(i).m)),1); %mad observed melt rate
    meltdata(1,12) = nanmedian(nanmean([melt(i).mHJ; melt(i).mWC],1)); %median parameterized melt rate
    meltdata(1,13) = mad(nanmean([melt(i).mHJ; melt(i).mWC],1),1); %mad parameterized melt rate
    meltdata(1,14) = nanmedian(melt(i).mHJ - melt(i).mWC); %median difference between HJ & WC melt rates
    if ~isempty(melt(i).oceanTavg)
        meltdata(1,15) = nanmedian(melt(i).oceanTavg - temp_weighted'); %median difference between observed & SOSE temps
        meltdata(1,16) = nanmedian(nanmean([melt(i).oceanUWC_est; melt(i).oceanUHJ_est],1)); %median parameterized water speed
        meltdata(1,17) = mad(nanmean([melt(i).oceanUWC_est(~isnan(melt(i).m)); melt(i).oceanUHJ_est(~isnan(melt(i).m))],1),1); %mad parameterized water speed
    else
        meltdata(1,15) = -9999;
        meltdata(1,16) = -9999;
        meltdata(1,17) = -9999;
    end
    dlmwrite(combined_filename,meltdata,'delimiter','\t','precision','%.2f','-append');
    clear meltdata;
    
    %report observed:modeled temp & inferred:modeled velocity
    if ~isempty(melt(i).oceanTavg)
        disp(['Observed - Modeled temperature = ',num2str(round(nanmedian(melt(i).oceanTavg - temp_weighted'),2)),' +/- ',num2str(round(mad(melt(i).oceanTavg - temp_weighted',1),2))]);
        disp(['Inferred:Modeled velocity = ',num2str(round(nanmedian(nanmean([melt(i).oceanUWC_est; melt(i).oceanUHJ_est])./vel_weighted),2)),' +/- ',num2str(round(mad(nanmean([melt(i).oceanUWC_est; melt(i).oceanUHJ_est])./vel_weighted,1),2))]);
        velratio(i) = nanmedian(nanmean([melt(i).oceanUWC_est; melt(i).oceanUHJ_est])./vel_weighted);
    else
        velratio(i) = NaN;
    end
    clear *weight*
    
    %replace negative Parameterized melt rates with 0s
    melt(i).mWC(melt(i).mWC<0)=NaN; melt(i).mHJ(melt(i).mHJ<0)=NaN;
    
    %add observed temps vs melt rate to plots
    figure(fig3); subplot(subpl);
    if ~isempty(melt(i).oceanTavg)
%         plT(3) = plot(melt(i).oceanTavg,365*melt(i).m,'s','markerfacecolor','none','color',[1 0.5 0.1],'linewidth',1.5); hold on;
%         plT(4) = plot(melt(i).oceanTbottom,365*melt(i).m,'^','markerfacecolor','none','color',[1 0.5 0.1],'linewidth',1.5); hold on;
        lat_area = melt(i).Asub-melt(i).Asurf; bas_area = melt(i).Asurf;
        plT(2) = plot(((melt(i).oceanTavg.*lat_area)+(melt(i).oceanTbottom.*bas_area))./(lat_area+bas_area),365*melt(i).m,'x','markerfacecolor','none','color',[1 0.5 0.1],'linewidth',1.5,'markersize',9); hold on;
        clear lat_area bas_area;
    end
    %format
    set(gca,'fontsize',20); grid on; 
    pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) pos(3) 0.07]);
    %standardize limits
    if mod(plot_loc(i),2) == 0
        if plot_loc(i) <= 4
            set(gca,'xlim',[-2 0.5],'ylim',[0 22],'ytick',[0:10:22]);
        else
            set(gca,'xlim',[-2 0.5],'ylim',[0 9],'ytick',[0:4:9]);
        end
    else
        set(gca,'xlim',[-2 0.5],'ylim',[0 70],'ytick',[0:30:70]);
    end
    text(min(get(gca,'xlim'))+0.6*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),0.85*max(get(gca,'ylim')),[char(plot_names(i)),' ',char(melt(i).dispname)],'fontsize',16);
    %add axes labels
    if plot_loc(i) == 14
        xlabel(['Water temperature (',char(176),'C)'],'fontsize',20);
    elseif plot_loc(i) == 15
        xlabel(['Water temperature (',char(176),'C)'],'fontsize',20); 
        ylblT = ylabel('Iceberg melt rate (m yr^{-1})','fontsize',20);
    elseif plot_loc(i) == 1
        set(gca,'xticklabel',[]); 
%         legT = legend(plT,'SOSE_{avg}','SOSE_{bottom}','obs._{avg}','obs._{bottom}'); 
        legT = legend(plT,'SOSE','observed'); 
        set(legT,'fontsize',16,'location','northoutside','NumColumns',2); set(gca,'position',[pos(1) pos(2) pos(3) 0.07]);
    else
        set(gca,'xticklabel',[]);
    end
    
    %plot observed vs parameterized melt rates
    figure(fig2); subplot(8,2,plot_loc(i));
    pl(1) = plot(melt(i).mHJ,365*melt(i).m,'^','color',[0 0 0],'markerfacecolor','none','linewidth',1.5); hold on;
    pl(2) = plot(melt(i).mWC,365*melt(i).m,'^','color',[0.5 0.5 0.5],'markerfacecolor','none','linewidth',1.5); hold on;
    disp(['Observed:HJ melt rate = ',num2str(round(nanmedian((365*melt(i).m)./melt(i).mHJ'),2)),' +/- ',num2str(round(mad((365*melt(i).m)./melt(i).mHJ',1),2))]);
    disp(['Observed:WC melt rate = ',num2str(round(nanmedian((365*melt(i).m)./melt(i).mWC'),2)),' +/- ',num2str(round(mad((365*melt(i).m)./melt(i).mWC',1),2))]);
    meltratio(i) = nanmedian((365*melt(i).m)./nanmean([melt(i).mWC(j); melt(i).mHJ(j)]));
    %format
    set(gca,'fontsize',20); grid on; 
    pos = get(gca,'position'); set(gca,'position',[pos(1) pos(2) pos(3) 0.07]);
    %standardize limits
    if mod(plot_loc(i),2) == 0
        if plot_loc(i) <= 4
            set(gca,'ylim',[0 22],'ytick',[0:10:22],'xlim',[0 3]);
        else
            set(gca,'ylim',[0 9],'ytick',[0:4:9],'xlim',[0 3]);
        end
    else
        set(gca,'ylim',[0 70],'ytick',[0:30:70],'xlim',[0 12],'xtick',[0:4:12]);
    end
    text(0.6*max(get(gca,'xlim')),0.85*max(get(gca,'ylim')),[char(plot_names(i)),' ',char(melt(i).dispname)],'fontsize',16);
    %add axes labels
    if plot_loc(i) == 14
        xlabel('Parameterized iceberg melt rate (m yr^{-1})','fontsize',20);
    elseif plot_loc(i) == 15
        xlabel('Parameterized iceberg melt rate (m yr^{-1})','fontsize',20); 
        ylbl = ylabel('Iceberg melt rate (m yr^{-1})','fontsize',20);
    elseif plot_loc(i) == 1
        set(gca,'xticklabel',[]); leg = legend(pl,'Holland & Jenkins','Weeks & Campbell'); set(leg,'fontsize',16);
    else
        set(gca,'xticklabel',[]);
    end
    
    %save the site-specific figure
    figure(site_fig); 
    subplot(1,3,1);
    legyr = legend(plyr,num2str([2013:1:2020]')); set(legyr,'position',[0.915 0.3 0.06 0.45]);
    set(gca,'ylim',[0 ceil(max([365*melt(i).m; melt(i).mHJ'; melt(i).mWC']))]);
    text(min(get(gca,'xlim'))+0.9*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.95*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'a)','fontsize',16);
    subplot(1,3,2);
    set(gca,'ylim',[0 ceil(max([365*melt(i).m; melt(i).mHJ'; melt(i).mWC']))]);
    rectangle('position',[min(get(gca,'xlim'))+0.025*(max(get(gca,'xlim'))-min(get(gca,'xlim'))) min(get(gca,'ylim'))+0.85*(max(get(gca,'ylim'))-min(get(gca,'ylim'))) 0.53*(max(get(gca,'xlim'))-min(get(gca,'xlim'))) 0.13*(max(get(gca,'ylim'))-min(get(gca,'ylim')))]);
    plot(min(get(gca,'xlim'))+0.05*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.948*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'^k','markerfacecolor','w','markersize',12); hold on;
    plot(min(get(gca,'xlim'))+0.05*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.89*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'xk','markerfacecolor','w','markersize',14,'linewidth',1.5); hold on;
    text(min(get(gca,'xlim'))+0.08*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.95*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'SOSE temperature','fontsize',16);
    text(min(get(gca,'xlim'))+0.08*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.89*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'observed temperature','fontsize',16);
    text(min(get(gca,'xlim'))+0.9*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.95*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'b)','fontsize',16);
    subplot(1,3,3);
    if ~isempty(melt(i).oceanTavg)
        set(gca,'ylim',[0 max([melt(i).oceanUHJ_est melt(i).oceanUWC_est(j)])+0.1*max([melt(i).oceanUHJ_est melt(i).oceanUWC_est(j)])]);
    else
        set(gca,'ylim',[0 max(get(gca,'ylim'))+0.1*max(get(gca,'ylim'))]);
    end
    rectangle('position',[min(get(gca,'xlim'))+0.025*(max(get(gca,'xlim'))-min(get(gca,'xlim'))) min(get(gca,'ylim'))+0.85*(max(get(gca,'ylim'))-min(get(gca,'ylim'))) 0.43*(max(get(gca,'xlim'))-min(get(gca,'xlim'))) 0.13*(max(get(gca,'ylim'))-min(get(gca,'ylim')))]);
    plot(min(get(gca,'xlim'))+0.05*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.948*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'^k','markerfacecolor','w','markersize',12); hold on;
    plot(min(get(gca,'xlim'))+0.05*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.89*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'xk','markerfacecolor','w','markersize',14,'linewidth',1.5); hold on;
    text(min(get(gca,'xlim'))+0.08*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.95*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'SOSE velocity','fontsize',16);
    text(min(get(gca,'xlim'))+0.08*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.89*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'inferred velocity','fontsize',16);
    text(min(get(gca,'xlim'))+0.9*(max(get(gca,'xlim'))-min(get(gca,'xlim'))),min(get(gca,'ylim'))+0.95*(max(get(gca,'ylim'))-min(get(gca,'ylim'))),'c)','fontsize',16);
    %save
    saveas(gcf,[char(melt(i).dispname),'-iceberg-melt-scatterplots.eps'],'epsc'); saveas(gcf,[char(melt(i).dispname),'-iceberg-melt-scatterplots.png'],'png');
    drawnow;
end
cd /Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt
% save('Antarctic-icebergmelt-comparison.mat','melt','-v7.3');

%format and save figures
cd /Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt
figure(fig1);
set(gca,'fontsize',20,'ytick',[1:360:length(SOSEx)],'yticklabel',[0:60:360],...
    'xtick',[1:48:length(SOSEy)],'xticklabel',[round(min(SOSEy)):3:round(max(SOSEy))]);
xlabel(['Latitude (',char(176),'N)'],'fontsize',20); ylabel(['Longitude (',char(176),'E)'],'fontsize',20); 
saveas(gcf,'Antarctic-SOSE-locations.eps','epsc'); saveas(gcf,'Antarctic-SOSE-locations.png','png');
figure(fig2);
set(ylbl,'position',[-2    350.0000   -1.0000]);
set(leg,'position',[0.61 0.09 0.2212 0.0446]);
saveas(gcf,'Antarctic-iceberg-melt-observed-vs-parameterized-scatterplots.eps','epsc'); saveas(gcf,'Antarctic-iceberg-melt-observed-vs-parameterized-scatterplots.png','png');
figure(fig3);
set(ylblT,'position',[-2.5    350.0000   -1.0000]);
set(legT,'position',[.59 0.08 0.2994 0.0603]);
saveas(gcf,'Antarctic-iceberg-melt-vs-temps-scatterplots.eps','epsc'); saveas(gcf,'Antarctic-iceberg-melt-vs-temps-scatterplots.png','png');
clear leg* pl plT ylbl*;


%% create geographically-arranged subplots that show the difference in observed & parameterized iceberg melt rates
close all; drawnow;

%set-up figure for subplots
melt_fig = figure; set(gcf,'position',[50 50 800 1000]); 
vel_fig = figure; set(gcf,'position',[450 50 800 1000]); 
yr_cmap = cmocean('dense',9); yr_cmap = yr_cmap(2:end,:);

%loop through melt structure & plot figures
for i = 1:length(melt)
    
    %identify the envelope for observed melt rate vs draft 
%     min_dmdd = min((365*melt(i).m)./melt(i).d); max_dmdd = max((365*melt(i).m)./melt(i).d); 
    DT = delaunayTriangulation(melt(i).d,(365*melt(i).m));
    C = convexHull(DT);
    
    %melt rate subplots
    figure(melt_fig);
    subpl = subplot(8,2,plot_loc(i)); set(gca,'fontsize',20); 
    plotpos = get(subpl,'position');
%     fill([0 max(melt(i).d) max(melt(i).d) 0], [0 max_dmdd*max(melt(i).d) min_dmdd*max(melt(i).d) 0],[0.5 0.5 0.5],'edgecolor','none','facealpha',0.5);
    fill(DT.Points(C,1),DT.Points(C,2),[0.5 0.5 0.5],'edgecolor','none','facealpha',0.5); hold on;
    if plot_loc(i) == 15
        for j = 1:length(yr_cmap)
            plyr(j) = plot(melt(i).d(1),nanmean([melt(i).mWC(1); melt(i).mHJ(1)]),'s','markerfacecolor',yr_cmap(j,:),'linewidth',1.5,'markeredgecolor','none');
        end
    end
    for j = 1:length(melt(i).d)
        plot(melt(i).d(j),nanmean([melt(i).mWC(j); melt(i).mHJ(j)]),'^k','markerfacecolor',yr_cmap(round(nanmean([melt(i).to(j) melt(i).tf(j)]))-2012,:),'markersize',8); hold on;
        if ~isempty(melt(i).oceanTavg)
            plot(melt(i).d(j),nanmean([melt(i).mWC_TSobs(j); melt(i).mHJ_TSobs(j)]),'x','color',yr_cmap(round(nanmean([melt(i).to(j) melt(i).tf(j)]))-2012,:),'markersize',9,'linewidth',1.5); hold on;
        end
    end
    %set axes
    grid on;
    if mod(plot_loc(i),2) == 0
        set(gca,'xlim',[0 400],'xticklabel',[],'fontsize',16);
        if plot_loc(i) >= 6
            set(gca,'ylim',[0 9],'ytick',[0:4:9]);
        else
            set(gca,'ylim',[0 22],'ytick',[0:10:22]);
        end
    else
        set(gca,'xlim',[0 800],'ylim',[0 70],'ytick',[0:30:70],'xticklabel',[],'fontsize',16);
    end
    %label bottom subplots
    if plot_loc(i) == 14
        set(gca,'xticklabel',[0:100:400]);
        xlabel('Iceberg draft (m b.s.l.)','fontsize',20);
    elseif plot_loc(i) == 15
        set(gca,'xticklabel',[0:200:800]);
        legyr = legend(plyr,num2str([2013:1:2020]'),'NumColumns',4);
        legpos = get(legyr,'position'); set(legyr,'position',[0.56 0.11 legpos(3) legpos(4)]); clear legpos;
        xlabel('Iceberg draft (m b.s.l.)','fontsize',20); 
        ylbl = ylabel('Parameterized iceberg melt rate (m yr^{-1})','fontsize',20);
        set(ylbl,'position',[-100 350 -1]);
    end
    %adjust subplot positions
    if mod(plot_loc(i),2) ~= 0
        set(gca,'position',[plotpos(1)-0.03 plotpos(2)+0.04/(plot_loc(i)+1) 1.1*plotpos(3) 1.1*plotpos(4)]);
    else
        set(gca,'position',[plotpos(1) plotpos(2)+0.04/(plot_loc(i)) 1.1*plotpos(3) 1.1*plotpos(4)]);
    end
    text(0.6*max(get(gca,'xlim')),0.85*max(get(gca,'ylim')),[char(plot_names(i)),' ',char(melt(i).dispname)],'fontsize',16);
    clear plotpos;
    
    %parameterized relative speed subplots
    figure(vel_fig); subpl2 = subplot(8,2,plot_loc(i));
    set(gca,'fontsize',20); grid on;
    for j = 1:length(melt(i).d)
        
        %identify the nearest SOSE cell with data
        [lon,lat] = ps2wgs(melt(i).x(j),melt(i).y(j),'StandardParallel',-71,'StandardMeridian',0);
        lon(lon<0) = lon(lon<0)+360;
        %         figure(fig1); plot(lat,lon,'xr','linewidth',2,'markersize',20); hold on;
        diffmap = sqrt((lon-SOSExgrid).^2 + (lat-SOSEygrid).^2);
        diffmap(Velmap==0 | isnan(Velmap)) = NaN;
        map_ind = find(diffmap == min(min(diffmap)));
        [x_ind,y_ind] = ind2sub(size(Tmap),map_ind);
        clear lon lat diffmap map_ind;
        
        %find the SOSE indices for the appropriate date range
        if melt(i).to(j) > 2013 && melt(i).tf(j) < 2020
            earlyref = find(abs(melt(i).to(j)-SOSEdate) == min(abs(melt(i).to(j)-SOSEdate)));
            lateref = find(abs(melt(i).tf(j)-SOSEdate) == min(abs(melt(i).tf(j)-SOSEdate)));
            Tprofs = squeeze(SOSET(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:lateref));
            Sprofs = squeeze(SOSES(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:lateref));
            Uprofs = squeeze(SOSEU(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:lateref));
            Vprofs = squeeze(SOSEV(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:lateref));
            
            clear earlyref lateref;
            
        elseif melt(i).to(j) < 2013
            %identify the partial years before 2013
            partyrs = ceil(12*(ceil(melt(i).to(j)) - melt(i).to(j)));
            
            %determine if the entire date range is before 2013
            if melt(i).tf(j) < 2013
                %identify the number of full years before 2013
                fullyrs = floor(melt(i).tf(j))-ceil(melt(i).to(j));
                
                %identify the end reference
                lateref = ceil(12*(melt(i).tf(j) - floor(melt(i).tf(j))));
                
                %combine averaged data for the full time period
                Tprofs = horzcat(squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                Sprofs = horzcat(squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                Uprofs = horzcat(squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                Vprofs = horzcat(squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                
            else
                %identify the number of full years before 2013
                fullyrs = 2013-ceil(melt(i).to(j));
                
                %identify the end reference
                lateref = find(abs(melt(i).tf(j)-SOSEdate) == min(abs(melt(i).tf(j)-SOSEdate)));
                
                %combine averaged data from before 2013 & post-2013 data
                Tprofs = horzcat(squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSET(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                Sprofs = horzcat(squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSES(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                Uprofs = horzcat(squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEU(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                Vprofs = horzcat(squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-partyrs:12)),...
                    fullyrs*squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEV(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:lateref)));
                
            end
            clear partyrs fullyrs lateref;
            
        elseif melt(i).tf(j) > 2019
            %identify the partial years after 2019
            partyrs = ceil(12*(melt(i).tf(j) - floor(melt(i).tf(j))));
            
            %determine if the entire date range is after 2019
            if melt(i).to(j) > 2019
                %identify the number of full years after 2019
                fullyrs = floor(melt(i).tf(j))-ceil(melt(i).to(j));
                
                %identify the start reference
                earlyref = ceil(12*(ceil(melt(i).to(j)) - melt(i).to(j)));
                
                %combine averaged data for the full time period
                Tprofs = horzcat(squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-earlyref:12)),...
                    fullyrs*squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                Sprofs = horzcat(squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-earlyref:12)),...
                    fullyrs*squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                Uprofs = horzcat(squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-earlyref:12)),...
                    fullyrs*squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                Vprofs = horzcat(squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),12-earlyref:12)),...
                    fullyrs*squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                
            else
                %identify the number of full years after 2019
                fullyrs = floor(melt(i).tf(j))-2019;
                
                %identify the start reference
                earlyref = find(abs(melt(i).to(j)-SOSEdate) == min(abs(melt(i).to(j)-SOSEdate)));
                
                %combine averaged data from before 2013 & post-2013 data
                Tprofs = horzcat(squeeze(SOSET(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:end)),...
                    fullyrs*squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSETmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                Sprofs = horzcat(squeeze(SOSES(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:end)),...
                    fullyrs*squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSESmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                Uprofs = horzcat(squeeze(SOSEU(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:end)),...
                    fullyrs*squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEUmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                Vprofs = horzcat(squeeze(SOSEV(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),earlyref:end)),...
                    fullyrs*squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),:)),...
                    squeeze(SOSEVmean(x_ind,y_ind,1:find(SOSEd<=melt(i).d(j),1,'last'),1:partyrs)));
                
            end
            clear partyrs fullyrs earlyref;
            
        end
        
        %average the profiles over time
        Tavg=nanmean(Tprofs,2); Savg=nanmean(Sprofs,2); Uavg=nanmean(Uprofs,2); Vavg=nanmean(Vprofs,2); davg = SOSEd(1:length(Tavg));
        Tavg(Savg==0) = []; Uavg(Savg==0) = []; Vavg(Savg==0) = []; davg(Savg==0) = []; Savg(Savg==0) = [];
        clear *profs;
        
        %determine the draft range that gives iceberg velocities that
        %minimize the difference between observed and modeled melt rates
        for k = 1:length(Tavg)
            depthref = k;
            
            %estimate iceberg velocity as the depth-averaged water velocity
            U_weight = nansum(Uavg(1:depthref).*SOSEH(1:depthref))./nansum(SOSEH(1:depthref));
            V_weight = nansum(Vavg(1:depthref).*SOSEH(1:depthref))./nansum(SOSEH(1:depthref));
            
            %estimate melt rate using the Holland & Jenkins/Silva melt equations
            [Mhj, T_sh, T_fp,a] = melt_forcedwater(Tavg,Savg,davg, sqrt((Uavg-U_weight).^2 + (Vavg-V_weight).^2), abs(melt(i).bergT));
            %calculate the average melt rate, weighting according to lateral & basal areas
            lat_area = melt(i).Asub(j)-melt(i).Asurf(j); bas_area = melt(i).Asurf(j);
            Mhj_weighted = nansum(Mhj.*SOSEH(1:length(Tavg)))./nansum(SOSEH(1:length(Tavg)));
            mHJ(k) = 31536000*((Mhj_weighted*lat_area)+(Mhj(end)*bas_area))./(lat_area+bas_area); %area-averaged melt rate (m yr^{-1})
            %estimate melt rate using the Weeks & Campbell melt equation
            Mwc = 0.037*((rho_sw/melt(i).rho(j))*diffM^(-7/15)*diffT^(2/3)*(c/L))*((sqrt((Uavg-U_weight).^2 + (Vavg-V_weight).^2).^(0.8).*(Tavg-T_fp))./((2*sqrt(melt(i).Asurf(j)/pi)).^0.2));
            Mwc_weighted = nansum(Mwc.*SOSEH(1:length(Tavg)))./nansum(SOSEH(1:length(Tavg)));
            mWC(k) = 31536000*((Mwc_weighted*lat_area)+(Mwc(end)*bas_area))./(lat_area+bas_area); %area-averaged melt rate (m yr^{-1})
            clear Mhj* Mwc* U_weight V_weight;
        end
        depthref = find(melt(i).m(j)./nanmean([mWC; mHJ]) == min(melt(i).m(j)./nanmean([mWC; mHJ])));
        clear mWC mHJ;
        %estimate iceberg velocity as the depth-averaged water velocity
        U_weighted(j) = nansum(Uavg(1:depthref).*SOSEH(1:depthref))./nansum(SOSEH(1:depthref));
        V_weighted(j) = nansum(Vavg(1:depthref).*SOSEH(1:depthref))./nansum(SOSEH(1:depthref));
        
        %plot the weighted velocities
        subplot(subpl2); plotpos2 = get(subpl2,'position');
        vel_weighted(j) = nansum(sqrt((Uavg-U_weighted(j)).^2 + (Vavg-V_weighted(j)).^2).*SOSEH(1:length(Tavg)))./nansum(SOSEH(1:length(Tavg)));
        if j == 1 && plot_loc(i) == 15
            for k = 1:length(yr_cmap)
                plyr2(k) = plot(melt(i).d(1),vel_weighted(1),'s','markerfacecolor',yr_cmap(k,:),'linewidth',1.5,'markeredgecolor','none'); hold on;
            end
        end
        plot(melt(i).d(j),vel_weighted(j),'^k','markerfacecolor',yr_cmap(round(nanmean([melt(i).to(j) melt(i).tf(j)]))-2012,:),'markersize',8); hold on;
        if ~isempty(melt(i).oceanTavg)
            plot(melt(i).d(j),nanmean([melt(i).oceanUHJ_est(j) melt(i).oceanUWC_est(j)]),'x','color',yr_cmap(round(nanmean([melt(i).to(j) melt(i).tf(j)]))-2012,:),'markersize',9,'linewidth',1.5); hold on;
        end
        
    end
    
    %set axes
    figure(vel_fig);  subplot(subpl2);
    grid on; set(gca,'ylim',[0 0.26],'fontsize',16);
    if mod(plot_loc(i),2) == 0
        set(gca,'xlim',[0 400],'xticklabel',[],'fontsize',16);
    else
        set(gca,'xlim',[0 800],'xticklabel',[],'fontsize',16);
    end
    %label bottom subplots
    if plot_loc(i) == 14
        set(gca,'xticklabel',[0:100:400],'ytick',[0:0.12:0.26]);
        xlabel('Iceberg draft (m b.s.l.)','fontsize',20);
    elseif plot_loc(i) == 15
        set(gca,'xticklabel',[0:200:800],'ylim',[0 2],'ytick',[0.9:2]);
        legyr2 = legend(plyr2,num2str([2013:1:2020]'),'NumColumns',4);
        legpos = get(legyr2,'position'); set(legyr2,'position',[0.56 0.11 legpos(3) legpos(4)]); clear legpos;
        xlabel('Iceberg draft (m b.s.l.)','fontsize',20); 
        ylb = ylabel('Relative velocity (m s^{-1})','fontsize',20);
        set(ylb,'position',[-100 4 -1]);
    elseif plot_loc(i) == 2 || plot_loc(i) == 8
        set(gca,'ylim',[0 2],'ytick',[0.9:2]);
    else
        set(gca,'ytick',[0:0.12:0.26]);
    end
    %adjust subplot positions
    if mod(plot_loc(i),2) ~= 0
        set(gca,'position',[plotpos2(1)-0.03 plotpos2(2)+0.04/(plot_loc(i)+1) 1.1*plotpos2(3) 1.1*plotpos2(4)]);
    else
        set(gca,'position',[plotpos2(1) plotpos2(2)+0.04/(plot_loc(i)) 1.1*plotpos2(3) 1.1*plotpos2(4)]);
    end
    text(0.6*max(get(gca,'xlim')),0.85*max(get(gca,'ylim')),[char(plot_names(i)),' ',char(melt(i).dispname)],'fontsize',16);
    clear plotpos; clear subpl*;
  
end

%save the figures
cd /Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt
figure(melt_fig);
saveas(gcf,'Antarctic-observed-and-parameterized-iceberg-melt-vs-draft-scatterplots.eps','epsc'); 
saveas(gcf,'Antarctic-observed-and-parameterized-iceberg-melt-vs-draft-scatterplots.png','png');
figure(vel_fig);
saveas(gcf,'Antarctic-inferred-water-velocity-vs-draft-scatterplots.eps','epsc'); 
saveas(gcf,'Antarctic-inferred-water-velocity-vs-draft-scatterplots.png','png');
disp('Saved gegraphically-arranged iceberg melt parameterization subplots');
