function [T] = plot_export_iceberg_melt_data(SL,dir_output,region_abbrev,DEM1,DEM2,plot_flag,table_flag)
% Function to plot iceberg melt data: meltwater flux vs. submerged area,
% melt rate vs. draft, and melt rate vs. sea level correction (i.e.,
% adjustment needed iceberg elevations so that they are with respect to a
% local sea level of zero meters)
% Ellyn Enderlin, Summer 2022
%
% INPUTS:   SL              structure containing meltwater flux, melt
%                           rate, and iceberg geometry data
%           DEM1            structure variable containing earlier DEM info
%           DEM2            structure variable containing later DEM info
%           iceberg_refs    number of iceberg to detect elevation change
%           dir_output      directory where all output files will be placed
%           region_abbrev   abbrevation of region used in file names
%           plot_flag       binary flag specifying plot generation (1=true)
%           table_flag      binary flag specifying table export (1=true)
%
% OUTPUTS:  If plot_flag = 1: Plots of meltwater flux vs. submerged area,
%           melt rate vs. draft, and melt rate vs. sea level correction.
%           If table_flag = 1: Summary data for iceberg location, size, 
%           meltwater flux, and melt rate exported to a csv table.


%plot
if plot_flag == 1
    %3-panel subplot: volume flux vs submerged area, melt rate vs draft,
    %melt rate vs sea level correction 
    dVdt = []; Asub = []; rho = []; H = []; m = []; coreg_zo = []; coreg_zf = []; berg_no =[]; SLref = [];
    for i = 1:length(SL)
        if SL(i).mean.dVdt > 0 && ~isempty(SL(i).mean.TA)
            dVdt = [dVdt SL(i).mean.dVdt];
            Asub = [Asub SL(i).mean.TA];
            rho = [rho (SL(i).initial.density+SL(i).final.density)/2];
            H = [H SL(i).mean.H];
            m = [m SL(i).mean.dHdt];
            coreg_zo = [coreg_zo SL(i).initial.coreg_z]; coreg_zf = [coreg_zf SL(i).final.coreg_z];
            berg_no = [berg_no; SL(i).name(end-1:end)];
            SLref = [SLref; i];
        end
    end
    figure; set(gcf,'position',[100 500 1500 600]);
    subplot(1,3,1);
    plot(Asub,dVdt,'ok','markersize',24,'markerfacecolor','w'); hold on;
    xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
    set(gca,'fontsize',20); xlabel('Submerged area (m^2)','fontsize',20); ylabel('Meltwater flux (m^3/d)','fontsize',20);
    for i = 1:length(m)
        text(double(Asub(i))-0.02*range(xlims),double(dVdt(i)),berg_no(i,:))
    end
    grid on;
    subplot(1,3,2);
    plot(H,m,'ok','markersize',24,'markerfacecolor','w'); hold on;
    xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
    set(gca,'fontsize',20); xlabel('Average iceberg thickness (m)','fontsize',20); ylabel('Melt rate (m/d)','fontsize',20);
    for i = 1:length(m)
        text(double(H(i))-0.02*range(xlims),double(m(i)),berg_no(i,:))
    end
    grid on;
    subplot(1,3,3);
    plot(coreg_zo-coreg_zf,dVdt,'ok','markersize',24,'markerfacecolor','w'); hold on;
    xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
    set(gca,'fontsize',20); xlabel('\Delta sea-level adjustment (m)','fontsize',20); ylabel('Meltwater flux (m^3/d)','fontsize',20);
    for i = 1:length(m)
        text(double(coreg_zo(i)-coreg_zf(i))-0.02*range(xlims),double(dVdt(i)),berg_no(i,:))
    end
    grid on;
    saveas(gcf,[dir_output,'/',DEM1.time,'-',DEM2.time,'/',region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt_scatterplots.eps'],'epsc');
    
    %"fix" melt rate estimates with bad local sea level adjustments
    sldz = coreg_zo-coreg_zf;  %sea level adjustment
    median_sldz = nanmedian(sldz); sldz_outlier = isoutlier(sldz);
    if sum(sldz_outlier) > 0
        disp('Automatically adjusting fluxes & melt rates for icebergs with clearly bad sea level estimates');
        median_slzo = nanmedian(coreg_zo); median_slzf = nanmedian(coreg_zf);
        rho_sw = 1026;  %kg m^-3
        for i = 1:length(coreg_zo)
            %         if (coreg_zo(i)-coreg_zf(i)) > (median_sldz + 1) || (coreg_zo(i)-coreg_zf(i)) < (median_sldz - 1)
            if sldz_outlier(i) == 1
                %confirm that you want to adjust sea level
                answer = questdlg(['Do you want to apply the median sea level adjustment to iceberg',berg_no(i,:),'?'],...
                    'Sea level adjustment','1) Yes!','2) No!','1) Yes!');
                switch answer
                    case '1) Yes!'
                        %adjust the volume change estimate so that it is calculated
                        %assuming the median sea level correction is accurate locally
                        rho_f = rho(i);
                        dZ_mean = (rho_sw/(rho_sw-rho_f))*(SL(i).mean.dz-((coreg_zo(i)-coreg_zf(i))-median_sldz));
                        dV_mean = SL(SLref(i)).mean.SA*(dZ_mean+SL(SLref(i)).SMB+SL(SLref(i)).creep_dz);
                        
                        %convert volume change to flux & melt rate
                        to = SL(SLref(i)).initial.time; tf = SL(SLref(i)).final.time;
                        dt = datenum(tf(1:12),'yyyymmddHHMM') - datenum(to(1:12),'yyyymmddHHMM');
                        dVdt_mean = dV_mean/dt;
                        dHdt_mean = dVdt_mean/SL(SLref(i)).mean.TA;
                        
                        %replace sea level adjustments & corrected volume flux (dVdt_mean) and melt rate (dHdt_mean) in structure
                        SL(SLref(i)).mean.dVdt = dVdt_mean; SL(SLref(i)).mean.dHdt = dHdt_mean;
                        SL(SLref(i)).initial.coreg_z = median_slzo; SL(SLref(i)).final.coreg_z = median_slzf;
                        
                end
                
                
            end
        end
    end
    clear sldz*;
    
    %2-panel subplot: volume flux vs submerged area, melt rate vs draft
    dVdt = []; Asub = []; H = []; m = []; coreg_zo = []; coreg_zf = []; berg_no =[]; SLref = [];
    for i = 1:length(SL)
        if SL(i).mean.dVdt > 0 && ~isempty(SL(i).mean.TA)
            dVdt = [dVdt SL(i).mean.dVdt];
            Asub = [Asub SL(i).mean.TA];
            H = [H SL(i).mean.H];
            m = [m SL(i).mean.dHdt];
            coreg_zo = [coreg_zo SL(i).initial.coreg_z]; coreg_zf = [coreg_zf SL(i).final.coreg_z];
            berg_no = [berg_no; SL(i).name(end-1:end)];
            SLref = [SLref; i]; 
        end
    end
    figure; set(gcf,'position',[100 100 1200 600]);
    subplot(1,2,1);
    plot(Asub,dVdt,'ok','markersize',24,'markerfacecolor','w'); hold on;
    xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
    set(gca,'fontsize',20); xlabel('Submerged area (m^2)','fontsize',20); ylabel('Meltwater flux (m^3/d)','fontsize',20);
    for i = 1:length(m)
        text(double(Asub(i))-0.02*range(xlims),double(dVdt(i)),berg_no(i,:))
    end
    grid on;
    subplot(1,2,2);
    plot(H,m,'ok','markersize',24,'markerfacecolor','w'); hold on;
    xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
    set(gca,'fontsize',20); xlabel('Average iceberg thickness (m)','fontsize',20); ylabel('Melt rate (m/d)','fontsize',20);
    for i = 1:length(m)
        text(double(H(i))-0.02*range(xlims),double(m(i)),berg_no(i,:))
    end
    grid on;
    disp('Iceberg meltwater flux should increase linearly with submerged area');
    disp('Iceberg melt rates should increase with thickness');
    clear xlims ylims;
    saveas(gcf,[dir_output,'/',DEM1.time,'-',DEM2.time,'/',region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt_scatterplots-adjusted.eps'],'epsc');
end

%save as a table in a text file
if table_flag == 1
    clear dt xo yo zo po Vo xf yf zf pf Vf coreg_z* dz* dVdt* draft* Asurf* Asub*;
    append_ref = 1;
    for i = 1:length(SL)
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
            disp(['Skipping over data for ',SL(i).name,' (',num2str(i),' in SL structure)']);
        end
    end
    bad_refs = find(draft<0); %remove data with negative thicknesses (unrealistic = error-prone)
    dt(bad_refs) = [];
    xo(bad_refs) = []; yo(bad_refs) = [];  zo(bad_refs) = []; Vo(bad_refs) = []; po(bad_refs) = []; coreg_zo(bad_refs) = [];
    xf(bad_refs) = []; yf(bad_refs) = [];  zf(bad_refs) = []; Vf(bad_refs) = []; pf(bad_refs) = []; coreg_zf(bad_refs) = [];
    dz(bad_refs) = []; dz_sigma(bad_refs) = []; dVdt(bad_refs) = []; dVdt_uncert(bad_refs) = [];
    draft(bad_refs) = []; draft_uncert(bad_refs) = []; Asurf(bad_refs) = []; Asurf_uncert(bad_refs) = []; Asub(bad_refs) = []; Asub_uncert(bad_refs) = [];
    clear bad_refs;
    T=table(dt',xo',yo',zo',po',Vo',xf',yf',zf',pf',Vf',coreg_zo',coreg_zf',dz',dz_sigma',dVdt',dVdt_uncert',draft',draft_uncert',Asurf',Asurf_uncert',Asub',Asub_uncert');
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
    writetable(T,[dir_output,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_meltinfo.csv']);
    disp('Text file written');
end


end