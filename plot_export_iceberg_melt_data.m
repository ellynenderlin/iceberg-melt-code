function [T] = plot_export_iceberg_melt_data(SL)

 %plot
    dVdt = []; Asub = []; H = []; m = []; coreg_zo = []; coreg_zf = []; berg_ref =[];
    for i = 1:length(SL)
        if SL(i).mean.dVdt > 0 && ~isempty(SL(i).mean.TA)
            dVdt = [dVdt SL(i).mean.dVdt];
            Asub = [Asub SL(i).mean.TA];
            H = [H SL(i).mean.H];
            m = [m SL(i).mean.dHdt];
            coreg_zo = [coreg_zo SL(i).initial.coreg_z]; coreg_zf = [coreg_zf SL(i).final.coreg_z];
            berg_ref = [berg_ref; SL(i).name(end-1:end)];
        end
    end
    figure; set(gcf,'position',[100 100 1200 600]);
    subplot(1,2,1);
    plot(Asub,dVdt,'ok','markersize',24,'markerfacecolor','w'); hold on;
    set(gca,'fontsize',20); xlabel('Submerged area (m^2)','fontsize',20); ylabel('Meltwater flux (m^3/d)','fontsize',20);
    for i = 1:length(m)
        text(double(Asub(i))-0.03e5,double(dVdt(i)),berg_ref(i,:))
    end
    grid on;
    subplot(1,2,2);
    plot(H,m,'ok','markersize',24,'markerfacecolor','w'); hold on;
    set(gca,'fontsize',20); xlabel('Average iceberg thickness (m)','fontsize',20); ylabel('Melt rate (m/d)','fontsize',20);
    for i = 1:length(m)
        text(double(H(i))-2,double(m(i)),berg_ref(i,:))
    end
    grid on;
    saveas(gcf,[dir_output,'/',DEM1.time,'-',DEM2.time,'/',region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt_scatterplots-adjusted.eps'],'epsc');
    
    %resave as a table in a text file
    cd(dir_iceberg);
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
%             SLref = strmatch(berg_ref(i,:),berg_nostring);
            disp(['Skipping over data for ',num2str(berg_ref(i,:))]);
        end
    end
    bad_refs = find(draft<0); %remove data with negative thicknesses (unrealistic = error-prone)
    dt(bad_refs) = []; 
    xo(bad_refs) = []; yo(bad_refs) = [];  zo(bad_refs) = []; Vo(bad_refs) = []; po(bad_refs) = []; coreg_zo(bad_refs) = []; 
    xf(bad_refs) = []; yf(bad_refs) = [];  zf(bad_refs) = []; Vf(bad_refs) = []; pf(bad_refs) = []; coreg_zf(bad_refs) = []; 
    dz(bad_refs) = []; dz_sigma(bad_refs) = []; dVdt(bad_refs) = []; dVdt_uncert(bad_refs) = []; 
    draft(bad_refs) = []; draft_uncert(bad_refs) = []; Asurf(bad_refs) = []; Asurf_uncert(bad_refs) = []; Asub(bad_refs) = []; Asub_uncert(bad_refs) = []; 
    clear bad_refs; 
    column_names = {'TimeSeparation (days)' 'X_i (m)' 'Y_i (m)' 'MedianZ_i (m)'...         'Density_i (kg/m^3)' 'Volume_i (m^3)' 'X_f (m)' 'Y_f (m)' 'MedianZ_f (m)'...         'Density_f (kg/m^3)' 'Volume_f (m^3)' 'VerticalAdjustment_i (m)' 'VerticalAdjustment_f (m)'...         'MeanElevationChange (m)' 'StdevElevationChange (m)' 'VolumeChangeRate (m^3/d)' 'VolumeChangeRate_uncert (m^3/d)'...         'MedianDraft_mean (m)' 'MedianDraft_range (m)' 'SurfaceArea_mean (m^2)' 'SurfaceArea_range (m^2)'...         'SubmergedArea_mean (m^3)','SubmergedArea_range (m^3)'};
        'Density_i (kg/m^3)' 'Volume_i (m^3)' 'X_f (m)' 'Y_f (m)' 'MedianZ_f (m)'...
        'Density_f (kg/m^3)' 'Volume_f (m^3)' 'VerticalAdjustment_i (m)' 'VerticalAdjustment_f (m)'...
        'MeanElevationChange (m)' 'StdevElevationChange (m)' 'VolumeChangeRate (m^3/d)' 'VolumeChangeRate_uncert (m^3/d)'...
        'MedianDraft_mean (m)' 'MedianDraft_range (m)' 'SurfaceArea_mean (m^2)' 'SurfaceArea_range (m^2)'...
        'SubmergedArea_mean (m^3)','SubmergedArea_uncert (m^3)'};
    T=table(dt',xo',yo',zo',po',Vo',xf',yf',zf',pf',Vf',coreg_zo',coreg_zf',dz',dz_sigma',dVdt',dVdt_uncert',draft',draft_uncert',Asurf',Asurf_uncert',Asub',Asub_uncert'); T.Properties.VariableNames = column_names;
    writetable(T,[dir_output,'/',DEM1.time,'-',DEM2.time,'/',region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_meltinfo.csv']);
    disp('Text file written');



end