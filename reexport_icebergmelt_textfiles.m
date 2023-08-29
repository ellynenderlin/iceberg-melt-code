%re-convert the iceberg melt mat files to tab-delimited text files with the full date range in the name
clearvars; close all; 
addpath('/Users/ellynenderlin/Research/miscellaneous/general-code/');
years = 2011:2023;

%navigate to the iceberg-melt directory
%Greenland (all files in one directory)
cd /Volumes/CALVING/Greenland_icebergs/iceberg_melt-backup/compiled-meltrates
region_name = 'Greenland'; region_abbrev = '';
%Antarctica (files in each region's directory)
% cd /Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt
% region_name = 'Widdowson-Biscoe'; region_abbrev = 'WG';
% disp(['Region: ',region_name]);
% cd(region_name);

%find all the files in the specified directory
datepairs = dir('*iceberg_melt.mat'); %if creating csvs for the first time
% datepairs = dir('*20*-20*.txt'); %if swapping csvs for txts
figure; set(gcf,'position',[100 100 1200 600]); sub1 = subplot(1,2,1); sub2 = subplot(1,2,2); drawnow;
iceberg_cmap = colormap(parula(length(years))); iceberg_cmap = iceberg_cmap(1:end-1,:);
subplot(sub1); for k = 1:length(iceberg_cmap); pl(k) = plot(-10000,-1,'+','markersize',20,'linewidth',2,'markeredgecolor',iceberg_cmap(k,:)); hold on; end
max_dVdt = []; max_Asub = []; max_melt = []; max_H = [];
for k = 1:length(datepairs)
    disp(datepairs(k).name);
    
    %automatically determine if loading a matfile or txtfile
    if contains(datepairs(k).name,'mat')
        load(datepairs(k).name); %if creating csvs for the first time
        append_ref = 1;
        for i = 1:length(SL)
            if isfield(SL,'mean') %format from 2020 to present
                if SL(i).mean.dVdt > 0 && ~isempty(SL(i).mean.TA)
                    dt(append_ref) = sum(SL(i).days);
                    xo(append_ref) = nanmean(SL(i).initial.x); yo(append_ref) = nanmean(SL(i).initial.y); zo(append_ref) = SL(i).initial.z_median; Vo(append_ref) = SL(i).initial.V;
                    xf(append_ref) = nanmean(SL(i).final.x); yf(append_ref) = nanmean(SL(i).final.y); zf(append_ref) = SL(i).final.z_median; Vf(append_ref) = SL(i).final.V;
                    po(append_ref)=SL(i).initial.density; pf(append_ref) = SL(i).final.density;
                    coreg_zo(append_ref) = SL(i).initial.coreg_z; coreg_zf(append_ref) = SL(i).final.coreg_z;
                    dz(append_ref) = SL(i).mean.dz; dz_sigma(append_ref) = SL(i).uncert.dz;
                    dVdt(append_ref) = SL(i).mean.dVdt; dVdt_uncert(append_ref) = max(SL(i).uncert.dVdt);
                    draft(append_ref) = SL(i).mean.draft; draft_uncert(append_ref) = SL(i).change.draft;
                    Asurf(append_ref) = SL(i).mean.SA; Asurf_uncert(append_ref) = SL(i).change.SA;
                    Asub(append_ref) = SL(i).mean.TA; Asub_uncert(append_ref) = SL(i).change.TA;
                    append_ref = append_ref+1;
                else
                    disp(['Skipping over data for ',num2str(append_ref),' (',num2str(i),' in SL structure)']);
                end
            else
                if SL(i).dVdt.mean > 0 && ~isempty(SL(i).TA.mean)
                    dt(append_ref) = sum(SL(i).days);
                    xo(append_ref) = nanmean(SL(i).xo); yo(append_ref) = nanmean(SL(i).yo); zo(append_ref) = median(SL(i).zo(~isnan(SL(i).zo))); Vo(append_ref) = SL(i).V.initial;
                    xf(append_ref) = nanmean(SL(i).xf); yf(append_ref) = nanmean(SL(i).yf); zf(append_ref) = median(SL(i).zf(~isnan(SL(i).zf))); Vf(append_ref) = SL(i).V.final;
                    po(append_ref)=900; pf(append_ref) = 900;
                    coreg_zo(append_ref) = SL(i).adjust_o; coreg_zf(append_ref) = SL(i).adjust_f;
                    if ~isfield(SL,'dz')
                        dz(append_ref) = SL(i).zo_median - SL(i).zf_median; dz_sigma(append_ref) = nanstd(SL(i).zo - SL(i).zf);
                        draft(append_ref) = 900/(1027-900).*nanmean([SL(i).zo_median,SL(i).zf_median]); 
                        draft_uncert(append_ref) = 900/(1027-900).*(SL(i).zo_median-SL(i).zf_median);
                    else
                        dz(append_ref) = SL(i).dz.mean; dz_sigma(append_ref) = SL(i).dz.stdev;
                        draft(append_ref) = SL(i).draft.mean; draft_uncert(append_ref) = SL(i).draft.initial - SL(i).draft.final;
                    end
                    dVdt(append_ref) = SL(i).dVdt.mean; dVdt_uncert(append_ref) = max(SL(i).dVdt.err);
                    Asurf(append_ref) = SL(i).SA.mean; Asurf_uncert(append_ref) = abs(SL(i).SA.initial - SL(i).SA.final)/2;
                    Asub(append_ref) = SL(i).TA.mean; Asub_uncert(append_ref) = abs(SL(i).TA.initial - SL(i).TA.final)/2;
                    append_ref = append_ref+1;
                else
                    disp(['Skipping over data for ',num2str(append_ref),' (',num2str(i),' in SL structure)']);
                end
            end
        end
        bad_refs = find(draft<0); %remove data with negative thicknesses (unrealistic = error-prone)
        dt(bad_refs) = [];
        xo(bad_refs) = []; yo(bad_refs) = [];  zo(bad_refs) = []; Vo(bad_refs) = []; po(bad_refs) = []; coreg_zo(bad_refs) = [];
        xf(bad_refs) = []; yf(bad_refs) = [];  zf(bad_refs) = []; Vf(bad_refs) = []; pf(bad_refs) = []; coreg_zf(bad_refs) = [];
        dz(bad_refs) = []; dz_sigma(bad_refs) = []; dVdt(bad_refs) = []; dVdt_uncert(bad_refs) = [];
        draft(bad_refs) = []; draft_uncert(bad_refs) = []; Asurf(bad_refs) = []; Asurf_uncert(bad_refs) = []; Asub(bad_refs) = []; Asub_uncert(bad_refs) = [];
        clear bad_refs;
    else
        M = dlmread(datepairs(k).name); %if swapping csvs for txts
        
        %filter out any bad data
        bad_refs = find(M(:,18)<0);
        M(bad_refs,:) = [];
        
        %extract data to plot
        %     dVdt = []; Asub = []; H = []; m = []; coreg_zo = []; coreg_zf = []; berg_ref =[];
        %     for i = 1:length(SL)
        %         if SL(i).mean.dVdt > 0 && ~isempty(SL(i).mean.TA)
        %             dVdt = [dVdt SL(i).mean.dVdt];
        %             Asub = [Asub SL(i).mean.TA];
        %             %extract the iceberg thickness
        %             if isfield(SL(i).mean,'H')
        %                 H = [H SL(i).mean.H];
        %             else
        %                 H = [H 1027/nanmean([SL(i).initial.density SL(i).final.density])*SL(i).mean.draft];
        %             end
        %             %extract the melt rate
        %             if isfield(SL(i).mean,'dHdt')
        %                 m = [m SL(i).mean.dHdt];
        %             else
        %                 m = [m SL(i).mean.dVdt./SL(i).mean.TA];
        %             end
        %         end
        %     end
        %     max_dVdt = [max_dVdt max(dVdt./86400)]; max_Asub = [max_Asub max(Asub)]; max_melt = [max_melt max(m)]; max_H = [max_H max(H)];
        dt = M(:,1);
        xo = M(:,2); yo = M(:,3); zo = M(:,4); po = M(:,5); Vo = M(:,6);
        xf = M(:,7); yf = M(:,8); zf = M(:,9); pf = M(:,10); Vf = M(:,11);
        coreg_zo = M(:,12); coreg_zf = M(:,13); dz = M(:,14); dz_sigma = M(:,15);
        dVdt = M(:,16); dVdt_uncert = M(:,17); draft = M(:,18); draft_uncert = M(:,19);
        Asurf = M(:,20); Asurf_uncert = M(:,21); Asub = M(:,22); Asub_uncert = M(:,23);
    end
    H = 1027./nanmean([po pf],2).*draft; m = dVdt./Asub;
    
    %extract the datestrings from the file name
    if contains(datepairs(k).name,region_name)
        d1 = datepairs(k).name(end-49:end-46); d2 = datepairs(k).name(end-34:end-31);
        output_name = [region_abbrev,datepairs(k).name(end-50:end-3),'csv'];
    else
        d1 = datepairs(k).name(4:11); d2 = datepairs(k).name(13:20); %assumes 2-letter abbreviation at file name beginning
        output_name = [datepairs(k).name(1:end-3),'info.csv'];
    end
    
    %plot the data
    subplot(sub1);
    plot(Asub,dVdt./86400,'+','markersize',20,'linewidth',2,'markeredgecolor',iceberg_cmap(round((convert_to_decimaldate(d1)+convert_to_decimaldate(d2))/2)-(years(1)-1),:)); hold on; 
    xlims = get(gca,'xlim'); set(gca,'xlim',[0 max(xlims)]);
    xticks = get(gca,'xtick'); set(gca,'xticklabel',xticks/10^6);
    set(gca,'fontsize',20); xlabel('Submerged area (km^2)','fontsize',20); ylabel('Meltwater flux (m^3/s)','fontsize',20);
    grid on;
    subplot(sub2);
    plot(H,m,'+','markersize',20,'linewidth',2,'markeredgecolor',iceberg_cmap(round((convert_to_decimaldate(d1)+convert_to_decimaldate(d2))/2)-(years(1)-1),:)); hold on;
    xlims = get(gca,'xlim'); set(gca,'xlim',[0 max(xlims)]);
    set(gca,'fontsize',20); xlabel('Average iceberg thickness (m)','fontsize',20); ylabel('Melt rate (m/d)','fontsize',20);
    grid on;
    drawnow;
    
    %resave as a tab-delimited text file
%     clear dt xo yo zo po Vo xf yf zf pf Vf coreg_z* dz* dVdt* draft* Asurf* Asub*;
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
    T=table(dt,xo,yo,zo,po,Vo,xf,yf,zf,pf,Vf,coreg_zo,coreg_zf,dz,dz_sigma,dVdt,dVdt_uncert,draft,2*draft_uncert,Asurf,2*Asurf_uncert,Asub,2*Asub_uncert);
    
    %save in the region directory
    writetable(T,output_name); clear T;
    
    clear SL M dt xo yo zo po Vo xf yf zf pf Vf coreg_zo coreg_zf dz dz_sigma dVdt dVdt_uncert draft draft_uncert Asurf Asurf_uncert Asub Asub_uncert;
    clear d1 d2 output_name;
end
subplot(sub1); leg = legend(pl,num2str([years]')); set(leg,'location','northwest');
ylims = get(gca,'ylim'); set(gca,'ylim',[0 max(ylims)]);
saveas(gcf,[region_name,'_iceberg-melt_subplots.png'],'png');