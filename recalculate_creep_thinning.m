%recalculate creep & adjust melt rate estimates as necessary
clear all; close all;
% region = 'LarsenA'; region_abbrev = 'LA';
% region = 'LarsenB'; region_abbrev = 'LB';
% region = 'Danco'; region_abbrev = 'DC';
% region = 'Palmer'; region_abbrev = 'PD';
% region = 'Biscoe'; region_abbrev = 'BS';
region = 'Marguerite'; region_abbrev = 'MA';
% region = 'Bugge'; region_abbrev = 'BI';

%cd to region directory & identify dates
cd /users/mariamadryak/Desktop/Antarctic_icebergs/
cd_to_region = ['cd ',region,'/']; eval(cd_to_region);
melt_dates = dir('20*-20*');
for i = 1:length(melt_dates)
    cd_to_date = ['cd ',melt_dates(i).name,'/iceberg_data/']; eval(cd_to_date);
    meltfile = dir('*_iceberg_melt.mat');
    disp(['Re-calculating creep for ',meltfile(1).name]);
    load_data = ['load ',meltfile(1).name,' SL']; eval(load_data);  
    
    %pull iceberg coordinates to re-extract RACMO air temps
    cd ../..
    meltinfo = dir('*iceberg_meltinfo.txt');
    lat = []; lon = []; psx = []; psy = [];
    for j = 1:length(meltinfo)
        mean_x = []; mean_y = []; mean_lat = []; mean_lon = [];
        plot_yrs = [str2num(meltinfo(j).name(end-37:end-34)) str2num(meltinfo(j).name(end-28:end-25))];
        load_data = ['M=dlmread(''',meltinfo(j).name,''');']; eval(load_data);
        xo = M(:,2); yo = M(:,3); %initial locations, median elev, density, volume
        xf = M(:,7); yf = M(:,8);  %same as above but final
        mean_x = nanmean([xo xf],2); mean_y = nanmean([yo yf],2);
        [mean_lon,mean_lat] = ps2wgs(mean_x,mean_y,'StandardParallel',-71,'StandardMeridian',0);
        lat=[lat; mean_lat]; lon=[lon; mean_lon]; psx = [psx; mean_x]; psy = [psy; mean_y];
        clear xo xf yo yf mean_*;
    end
    avg_lat = nanmean(lat); avg_lon = nanmean(lon); 
    %pull RACMO temps
    cd ../RACMO_Antarctica
    prompt = 'Are you looking at icebergs along the Antarctic Peninsula (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        AP=1;
        runoff_lat = ncread('RACMO2.3p2_XPEN055_runoff_daily_2011_2016.nc','lat');
        runoff_lon = ncread('RACMO2.3p2_XPEN055_runoff_daily_2011_2016.nc','lon');
        icetemp = ncread('RACMO2.3p2_XPEN055_T2m_monthly_1979_2016.nc','t2m'); icetemp(icetemp==-9999)=NaN; %assume the ice temp matches the avg long-term temp MINUS 10C (for creep estimates)
    else
        AP=0;
        runoff_lat = ncread('RACMO2.3p2_ANT27_runoff_daily_2011_2016.nc','lat');
        runoff_lon = ncread('RACMO2.3p2_ANT27_runoff_daily_2011_2016.nc','lon');
        icetemp = ncread('RACMO2.3p2_ANT27_T2m_monthly_1979_2016.nc','t2m'); icetemp(icetemp==-9999)=NaN; %assume the ice temp matches the avg long-term temp MINUS 10C (for creep estimates)
    end
    lat_diff = abs(avg_lat(1)*ones(size(runoff_lat)) - runoff_lat);
    lon_diff = abs(avg_lon(1)*ones(size(runoff_lon)) - runoff_lon);
    diff_map = sqrt(lat_diff.^2+lon_diff.^2);
    RACMO_ref = find(diff_map==min(min(diff_map)));
    [RACMOy RACMOx] = ind2sub(size(squeeze(nanmean(icetemp(:,:,1,270:450),4))),RACMO_ref);
    disp(['RACMO x-reference = ',num2str(RACMOx),' & y-reference = ',num2str(RACMOy)]);
    %adjust coordinates?
    disp('Adjust coordinates (if necessary) to extract surface melt estimates');
    figure; set(gcf,'position',[500 100 900 900]);
    runoff_cmap = colormap(jet(10001)); runoff_cmap(1,:) = [1 1 1]; colorbar;
    imagesc(nanmean(icetemp,4)-273.15); colormap(gca,runoff_cmap); hold on; set(gca,'ydir','reverse'); %plot coordinates on average melt from peak summer melt rates in the big 2012 melt season
    disp(['ice temperature at RACMO pixel = ',num2str(nanmean(icetemp(RACMOy,RACMOx,:,:))-273.15),' C']);
    cbar = colorbar; set(get(cbar,'ylabel'),'string','Average air temp. (C)');
    plot(RACMOx,RACMOy,'ok','markerfacecolor','k','markeredgecolor','w','markersize',20); hold on;
    plot(RACMOx,RACMOy,'xw','markersize',20); hold on;
    set(gca,'xlim',[RACMOx-10 RACMOx+10],'ylim',[RACMOy-10 RACMOy+10]);
    prompt = 'Modify coordinates if the marker is in a region with temps >0C. Do the coordinates need to be modified (y/n)?';
    str = input(prompt,'s');
    if strmatch(str,'y')==1
        disp('To change coordinates, identify new coordinates using tick marks and type RACMOx=XX; RACMOy=YY; dbcont ');
        keyboard
    end
    close(gcf); drawnow;
    for j = 1:length(SL)
        if length(SL(j).days) > 0
            SL(j).old_creep.airtemp = SL(j).airtemp; SL(j).old_creep.ratefactor = SL(j).ratefactor; %rename old air temp & rate factor
            if AP ==1
                SL(j).airtemp = nanmean(icetemp(RACMOy,RACMOx,:))-5;
            else
                SL(j).airtemp = nanmean(icetemp(RACMOy,RACMOx,:))-10; %pull new air temp to make sure it's correct
            end
            data_flag(j) = 1;
        else
            data_flag(j) = 0;
        end
    end
    clear RACMO* runoff* *lat* *lon* *diff* icetemp;
    sl = SL(data_flag==1); SL = sl; clear sl;
    
    %calculate precise time separation
    to = SL(1).initial.time; tf = SL(1).final.time;
    if mod(str2num(to(1:4)),4)==0; doyso=366; modayso = [31 29 31 30 31 30 31 31 30 31 30 31]; else doyso=365; modayso = [31 28 31 30 31 30 31 31 30 31 30 31]; end
    if mod(str2num(tf(1:4)),4)==0; doysf=366; modaysf = [31 29 31 30 31 30 31 31 30 31 30 31]; else doysf=365; modaysf = [31 28 31 30 31 30 31 31 30 31 30 31]; end
    doyo = sum(modayso(1:str2num(to(5:6))))-31+str2num(to(7:8)); doyf = sum(modaysf(1:str2num(tf(5:6))))-31+str2num(tf(7:8));
    if str2num(tf(1:4)) == str2num(to(1:4))
        ddays = doyf-doyo+1;
    elseif str2num(tf(1:4)) - str2num(to(1:4)) == 1
        ddays = doyf + (doyso-doyo)+1;
    else
        years = str2num(to(1:4)):1:str2num(tf(1:4));
        for k = 1:length(years)
            if mod(years(k),4)==0
                doys(k)=366;
            else
                doys(k) = 365;
            end
        end
        
        %calculate the sum of days differently if during a leap year
        if doyo > sum(modayso(1:2))
            ddays = doyf + sum(doys(2:end-1)) + (365-doyo)+1;
        else
            ddays = doyf + sum(doys(2:end-1)) + (366-doyo)+1;
        end
    end
    hrs_o = ((str2num(to(13:14))/(60*60*24))+(str2num(to(11:12))/(60*24))+(str2num(to(9:10))/24));
    hrs_f = ((str2num(tf(13:14))/(60*60*24))+(str2num(tf(11:12))/(60*24))+(str2num(tf(9:10))/24));
    dhrs = hrs_f - hrs_o;
    dt = ddays + dhrs; %DECIMAL DAYS
    
%     %find which ones included an adjusted (average) sea level correction
%     coreg_zo = []; coreg_zf = []; rho_sw = 1026;
%     for j = 1:length(SL)
%         if ~isempty(SL(j).mean.dHdt)
%             coreg_zo = [coreg_zo SL(j).initial.coreg_z]; coreg_zf = [coreg_zf SL(j).final.coreg_z];
%         end
%     end
%     median_sldz = nanmedian(coreg_zo-coreg_zf); median_slzo = nanmedian(coreg_zo); median_slzf = nanmedian(coreg_zf);
    
    
    %back-out volume change
    figure;
    for j = 1:length(SL)
        if ~isempty(SL(j).mean.dHdt) %icebergs already flagged as bad have empty dHdt fields
        dV(j) = SL(j).mean.dVdt.*dt;
        dH(j) = dV(j)./SL(j).mean.SA - SL(j).creep_dz; %same as (rho_sw/(rho_sw-nanmean([SL(j).initial.density SL(j).final.density])))*SL(j).mean.dz+SL(j).SMB
        
        %recalculate creep
        rho_f = nanmean([SL(j).initial.density SL(j).final.density]);
        rf_o = 3.5e-25; %rate factor at -10C (Pa^-3 s^-1) from Cuffey & Paterson p. 74
        if SL(j).airtemp <= 263
            rf=rf_o*exp((-60000/8.314)*((1/SL(j).airtemp)-(1/263))); %rate factor for cold ice
        else
            rf=rf_o*exp((-134000/8.314)*((1/SL(j).airtemp)-(1/263))); %rate factor for nearly-temperate ice
        end
        SL(j).ratefactor = rf; B = rf^(-1/3);
        creep(j) = (SL(j).mean.H*(-2*3^(-2)*((rho_f*9.81*SL(j).mean.z)/(2*B))^3)*(dt*86400)); %from Thomas et al. (1973) eqn 19 & assuming incompressibility
        disp(['Old creep vs new creep = ',num2str(SL(j).creep_dz),' vs ',num2str(creep(j))]);
        
        %rename/resave estimates with old creep values
        SL(j).old_creep.creep_dz = SL(j).creep_dz;
        SL(j).old_creep.range.dVdt(1) = SL(j).range.dVdt(1); SL(j).old_creep.range.dVdt(2) = SL(j).range.dVdt(2);
        SL(j).old_creep.uncert.dVdt = SL(j).uncert.dVdt;
        SL(j).old_creep.mean.dVdt = SL(j).mean.dVdt;
        SL(j).old_creep.range.dHdt(1) = SL(j).range.dHdt(1); SL(j).old_creep.range.dHdt(2) = SL(j).range.dHdt(2);
        SL(j).old_creep.uncert.dHdt = SL(j).uncert.dHdt;
        SL(j).old_creep.mean.dHdt = SL(j).mean.dHdt;
        
        %recalculate meltwater fluxes & melt rates
        SL(j).creep_dz = creep(j);
        %volume flux errors/uncertainties
        SL(j).mean.dVdt = ((dH(j) + creep(j))*SL(j).mean.SA)/dt; %m^3 d^-1 of iceberg-equiv flux
        dVdt_stdev = (((SL(j).old_creep.range.dVdt(2)-SL(j).old_creep.range.dVdt(1))/6)/SL(j).old_creep.mean.dVdt)*SL(j).mean.dVdt;
        dVdt_total_err = SL(j).old_creep.uncert.dVdt;
        SL(j).range.dVdt(1) = SL(j).mean.dVdt-3*dVdt_stdev; SL(j).range.dVdt(2) =SL(j).mean.dVdt+3*dVdt_stdev;
        SL(j).uncert.dVdt = (dVdt_total_err/SL(j).old_creep.mean.dVdt)*SL(j).mean.dVdt;
        %melt rate errors/uncertainties
        dHdt_stdev = (((SL(j).old_creep.range.dHdt(2)-SL(j).old_creep.range.dHdt(1))/6)/SL(j).old_creep.mean.dHdt)*SL(j).mean.dHdt;
        dHdt_total_err = SL(j).old_creep.uncert.dHdt;
        SL(j).mean.dHdt = SL(j).mean.dVdt/SL(j).mean.TA;
        SL(j).range.dHdt(1) = SL(j).mean.dHdt-3*dHdt_stdev; SL(j).range.dHdt(2) = SL(j).mean.dHdt+3*dHdt_stdev;
        SL(j).uncert.dHdt = (dHdt_total_err/SL(j).old_creep.mean.dHdt)*SL(j).mean.dHdt;
        
%         %if melt rate is negative, check (and adjust if appropriate) sea level correction
%         if SL(j).mean.dVdt < 0 && ((coreg_zo(j)-coreg_zf(j)) > (median_sldz + 1) || (coreg_zo(j)-coreg_zf(j)) < (median_sldz - 1))
%             dVdt_correction = (SL(j).mean.SA*((coreg_zo(j)-coreg_zf(j))-median_sldz))/dt;
%             
%             
%         end
        
        %plot the updated and old data
        plot(SL(j).mean.TA,SL(j).old_creep.mean.dVdt./86400,'ok','markerfacecolor','w'); hold on;
        plot(SL(j).mean.TA,SL(j).mean.dVdt./86400,'xr'); hold on;
        clear d*dt_stdev d*dt_total_err;
        drawnow;
        
        %remove bad bergs with negative melt rates
        if SL(j).mean.dVdt < 0
            SL(j).mean.TA = []; SL(j).mean.dHdt = [];
        end
        
        
        end
    end
    
    %resave 
    disp('saving melt rates');
    cd_to_date = ['cd ../',region,'/',melt_dates(i).name,'/iceberg_data/']; eval(cd_to_date);
    save_file = ['save(''',meltfile(1).name,''',','''SL''',',','''-v7.3''',')']; eval(save_file);
    clear dt xo yo zo po Vo xf yf zf pf Vf coreg_z* dz* dVdt* draft* Asurf* Asub*;
    column_names = ['dt','xo','yo','zo','po','Vo','xf','yf','zf','pf','Vf','coreg_zo','coreg_zf','dz','dz_sigma','dVdt','dVdt_uncert','draft','draft_uncert','Asurf','Asurf_uncert','Asub','Asub_uncert'];
    append_ref = 1;
    for j = 1:length(SL)
        if SL(j).mean.dVdt > 0 && ~isempty(SL(j).mean.TA)
            dt(append_ref) = sum(SL(j).days);
            xo(append_ref) = nanmean(SL(j).initial.x); yo(append_ref) = nanmean(SL(j).initial.y); zo(append_ref) = SL(j).initial.z_median; Vo(append_ref) = SL(j).initial.V;
            xf(append_ref) = nanmean(SL(j).final.x); yf(append_ref) = nanmean(SL(j).final.y); zf(append_ref) = SL(j).final.z_median; Vf(append_ref) = SL(j).final.V;
            po(append_ref)=SL(j).initial.density; pf(append_ref) = SL(j).final.density;
            coreg_zo(append_ref) = SL(j).initial.coreg_z; coreg_zf(append_ref) = SL(j).final.coreg_z;
            dz(append_ref) = SL(j).mean.dz; dz_sigma(append_ref) = SL(j).uncert.dz;
            dVdt(append_ref) = SL(j).mean.dVdt; dVdt_uncert(append_ref) = max(SL(j).uncert.dVdt);
            draft(append_ref) = SL(j).mean.draft; draft_uncert(append_ref) = SL(j).change.draft/2;
            Asurf(append_ref) = SL(j).mean.SA; Asurf_uncert(append_ref) = SL(j).change.SA/2;
            Asub(append_ref) = SL(j).mean.TA; Asub_uncert(append_ref) = SL(j).change.TA/2;
            append_ref = append_ref+1;
        else
            disp(['Skipping over data for ',num2str(j)]);
        end
    end
    column_vals = [dt' xo' yo' zo' po' Vo' xf' yf' zf' pf' Vf' coreg_zo' coreg_zf' dz' dz_sigma' dVdt' dVdt_uncert' draft' draft_uncert' Asurf' Asurf_uncert' Asub' Asub_uncert'];
    cd ../..
    dlmwrite([region_abbrev,'_',to(1:8),'-',tf(1:8),'_iceberg_meltinfo.txt'],column_vals,'delimiter','\t');
    disp('...saved!');

    
    clear SL rho_f B to doyo ddays years doys hrs_o hrs_f dhrs rf;
end
disp(['Done adjusting data for ',region]);


