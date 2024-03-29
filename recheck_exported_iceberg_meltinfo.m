%Check dVdt adjustments in end-product textfiles. Used to check:
%(1) RACMO SMB (account for swapped X&Y coords) 
%   - fixed in extract_RACMO_params.m
%(2) creep thinning calculations (account for 365x over-estimation of creep
%due to incorrect unit conversion) 
%   - fixed in convert_Antarctic_iceberg_elev_change_to_meltrates.m


%% initialize
clearvars; close all;
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code','/Users/ellynenderlin/Research/miscellaneous/general-code/cmocean');

%specify file locations
root_dir = '/Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt/';
RACMO_path = '/Users/ellynenderlin/Research/miscellaneous/RACMO2.3_Antarctica/';
dir_code = '/Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt-code/';
addpath(dir_code);

%specify study site names
region = [{'Edgeworth-LarsenA'},{'Crane-LarsenB'},{'Ronne'},{'Filchner'},{'Amery'},{'Totten'},{'Mertz'},...
    {'Thwaites'},{'Ferrigno-Eltanin'},{'Seller-Bugge'},{'Heim-Marguerite'},{'Widdowson-Biscoe'},{'Cadman-Palmer'},{'Blanchard-Danco'},{'Leonardo-Danco'}];
leg_names = [{'Edgeworth'},{'Crane'},{'Ronne'},{'Filchner'},{'Polar Times'},{'Totten'},{'Mertz'},{'Thwaites'},{'Ferrigno'},{'Seller'},{'Heim'},{'Widdowson'},{'Cadman'},{'Blanchard'},{'Leonardo'}];

%set densities: assume steady-state profiles (like the Ligtenberg et al. (2014) paper
%describing these data), as supported by fairly constant AIS climate over the past ~40 years
rho_sw = 1026;  rho_sw_err = 2; %kg m^-3
rho_i = 917; rho_i_err = 10; %kg m^-3
rf_o = 3.5e-25; %rate factor at -10C (Pa^-3 s^-1) from Cuffey & Paterson p. 74

%% check meltwater flux & melt rate estimates
warning off;
answer = questdlg('Where are you working?',...
    'Iceberg Location','1) Greenland','2) Antarctica','1) Greenland');
switch answer
    case '1) Greenland'
        geography = 0;
    case '2) Antarctica'
        geography = 1;
end

% loop through csvs for each site & pull RACMO SMB from coords
for j = 1:length(region)
    %navigate to the region folder
    cd([root_dir,char(region(j))]);
    csvs = dir('*meltinfo.csv');
    
    %identify the region where you are working for RACMO data extraction
    fprintf(['Site: ',char(region(j)),' \n']);
    
    %loop through each date & adjust volume change estimates
    for k = 1:length(csvs)
        cd([root_dir,char(region(j))]);
        T = readtable(csvs(k).name);
        Tp1 = table2array(T(:,1:15)); Tp2 = table2array(T(:,18:23));
        DEM1.time = csvs(k).name(4:17); DEM2.time = csvs(k).name(19:32);
        disp(['Date range = [',DEM1.time,',',DEM2.time,']']);
        %fix dates for csvs that contain data from multiple DEMs for the same date range
        if contains(DEM1.time,'X'); DEM1.time = [DEM1.time(1:end-1),'0']; end
        if contains(DEM2.time,'X'); DEM2.time = [DEM2.time(1:end-1),'0']; end
        
        %pull out needed variables for RACMO data extraction
        berg_x = nanmean([table2array(T(:,2)); table2array(T(:,7))]);
        berg_y = nanmean([table2array(T(:,3)); table2array(T(:,8))]);
        berg_dates = [DEM1.time; DEM2.time];
%         dt = datenum(DEM2.time,'yyyymmddHHMMSS') - datenum(DEM1.time,'yyyymmddHHMMSS');
        
        %extract air temp & firn density info from RACMO
        [dt,iceberg_avgtemp,surfmelt,firnair,density,f,ci] = extract_RACMO_params(RACMO_path,geography,berg_x,berg_y,berg_dates);
        fprintf('Surface melt = %4.3f m \n',surfmelt);
        fprintf('Time separation = %3.1f days \n',dt);
        
        %loop through each iceberg
        for l = 1:height(T)
            fprintf('Average density = %4.1f \n',nanmean([table2array(T(l,5)) table2array(T(l,10))]));
            %convert elevation change to thickness change
            dZ_mean = (rho_sw/(rho_sw-nanmean([table2array(T(l,5)) table2array(T(l,10))])))*table2array(T(l,14));
            dZ_stdev = (rho_sw/(rho_sw-nanmean([table2array(T(l,5)) table2array(T(l,10))])))*table2array(T(l,15));
            
            %correct for SMB loss
            dH_SMBadjust_mean = dZ_mean-surfmelt;
            
            %correct for creep thinning
            if iceberg_avgtemp <= 263
                rf=rf_o*exp((-60000/8.314)*((1/iceberg_avgtemp)-(1/263))); %rate factor for cold ice
            else
                rf=rf_o*exp((-134000/8.314)*((1/iceberg_avgtemp)-(1/263))); %rate factor for nearly-temperate ice
            end
            B = rf^(-1/3); %Pa s^1/3
            creep = ((-1/(2*sqrt(3)))*((nanmean([table2array(T(l,5)) table2array(T(l,10))])*9.81*nanmean([table2array(T(l,4)) table2array(T(l,9))]))/(2*sqrt(3)))^3*(1-(nanmean([table2array(T(l,5)) table2array(T(l,10))])/rho_sw))^3)/(B^3); %creep thinning rate (1/s)
            creep_dz = ((rho_sw/(rho_sw-nanmean([table2array(T(l,5)) table2array(T(l,10))]))).*nanmean([table2array(T(l,4)) table2array(T(l,9))]).*(86400*creep)*dt);
            fprintf('Elevation change w/o creep = %4.3f m \n',dH_SMBadjust_mean);
            fprintf('Creep thinning = %4.3f m \n',creep_dz);
            dH_submelt = dH_SMBadjust_mean + creep_dz; %integrate creep over the ice thickness & over the time period
            
            %recalculate volume change due to submarine melting
            dV_mean = table2array(T(l,20))*dH_submelt;
            
            %recalculate volume change rate due to submarine melting
            dVdt_mean = dV_mean./dt;
            dVdt_total_err = (table2array(T(l,17))./table2array(T(l,16))).*dVdt_mean;
            
            %display estimate of iceberg submarine melt volume flux before adjustment for SMB
            fprintf('Submarine meltwater flux (m^3/s): ');
            fprintf('Previous = %4.3f & Updated = %4.3f \n',table2array(T(l,16))./86400,dVdt_mean./86400); %volume change rate (m^3/d)
            fprintf('Submarine melt rate (m/yr): ');
            fprintf('Previous = %3.1f & Updated = %3.1f \n',365*table2array(T(l,16))./table2array(T(l,22)),365*dVdt_mean./table2array(T(l,22))); %volume change rate (m^3/d)
            
            %replace data dVdt & uncertainty
            dVdt(l,1) = dVdt_mean; dVdt_err(l,1) = dVdt_total_err;
            
            %clear variables for each iceberg
            clear d*_mean d*_stdev rf B creep* d*melt d*_total_err;
        end
        disp('close figure to advance');
        fig = figure; 
        if size(dVdt,1) == size(table2array(T(:,22)),1)
            plot(365*table2array(T(:,16))./table2array(T(:,22)),365*dVdt./table2array(T(:,22)),'xk','markersize',12); hold on;
        else
            plot(365*table2array(T(:,16))./table2array(T(:,22)),365*dVdt'./table2array(T(:,22)),'xk','markersize',12); hold on;
        end
        grid on;
        xlabel('Old melt rate (m/yr)'); ylabel('New melt rate (m/yr)');
        waitfor(fig);
        
        %remake the table
        clear T;
        T = table(repmat(dt,size(Tp1(:,2))),Tp1(:,2),Tp1(:,3),Tp1(:,4),Tp1(:,5),Tp1(:,6),Tp1(:,7),...
            Tp1(:,8),Tp1(:,9),Tp1(:,10),Tp1(:,11),Tp1(:,12),Tp1(:,13),Tp1(:,14),Tp1(:,15),...
            dVdt,dVdt_err,Tp2(:,1),Tp2(:,2),Tp2(:,3),Tp2(:,4),Tp2(:,5),Tp2(:,6));
        column_names = ["TimeSeparation", "X_i", "Y_i", "MedianZ_i",...
            "Density_i", "Volume_i", "X_f", "Y_f", "MedianZ_f",...
            "Density_f", "Volume_f", "VerticalAdjustment_i", "VerticalAdjustment_f",...
            "ElevationChange_mean", "ElevationChange_stdev", "VolumeChangeRate", "VolumeChangeRate_uncert",...
            "MedianDraft_mean", "MedianDraft_range", "SurfaceArea_mean", "SurfaceArea_range",...
            "SubmergedArea_mean","SubmergedArea_uncert"];
        column_units = ["days", "m", "m", "m",...
            "kg/m^3", "m^3", "m", "m", "m",...
            "kg/m^3", "m^3", "m", "m,"...
            "m", "m", "m^3/d", "m^3/d",...
            "m", "m", "m^2", "m^2",...
            "m^2", "m^2"];
        T.Properties.VariableNames = column_names; T.Properties.VariableUnits = column_units;
        writetable(T,[root_dir,char(region(j)),'/',csvs(k).name]);
        disp('Re-exported csv!');
        disp('moving on to next date ...');
        
        %clear variables
        clear iceberg_avgtemp surfmelt firnair density f ci;
        clear T* DEM1 DEM2 berg_* dVdt dVdt_err dt*;
    end
    clear csvs;
    disp('moving on to next site ...');
    disp(' ');
end




