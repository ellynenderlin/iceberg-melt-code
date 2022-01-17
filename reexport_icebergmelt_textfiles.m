%re-convert the iceberg melt mat files to tab-delimited text files with the full date range in the name
clearvars; close all; 

%navigate to the iceberg-melt directory
% cd /volumes/EllynEnderlin/common/Antarctica-icebergs/iceberg-melt
cd /Users/ellynenderlin/Research/NSF_Antarctic-Icebergs/iceberg-melt

%modify the region name to repeat this process for each location
region_name = 'Amery'; disp(['Region: ',region_name]);
cd_to_region = ['cd ',region_name]; eval(cd_to_region);

%find all the date directories
% datepairs = dir('20*-20*');
datepairs = dir('*20*-20*.txt');
figure; set(gcf,'position',[100 100 1200 600]); sub1 = subplot(1,2,1); sub2 = subplot(1,2,2);
iceberg_cmap = colormap(parula(9)); iceberg_cmap = iceberg_cmap(1:8,:);
subplot(sub1); for k = 1:length(iceberg_cmap); pl(k) = plot(-10000,-1,'+','markersize',20,'linewidth',2,'markeredgecolor',iceberg_cmap(k,:)); hold on; end
max_dVdt = []; max_Asub = []; max_melt = []; max_H = [];
for k = 1:length(datepairs)
    disp(datepairs(k).name);
%     cd_to_data = ['cd ',datepairs(k).name,'/iceberg_data/']; eval(cd_to_data);
%     matfile = dir('*iceberg_melt.mat'); load(matfile(1).name);
    M = dlmread(datepairs(k).name);
        
    %filter out any bad data
    bad_refs = find(M(:,18)<0);
    M(bad_refs,:) = [];
    dlmwrite(datepairs(k).name,M,'delimiter','\t');
    
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
    H = 1027./nanmean([po pf],2).*draft; m = dVdt./Asub;
    
    %plot the data
    subplot(sub1);
%     plot(Asub,dVdt./86400,'+','markersize',20,'linewidth',2,'markeredgecolor',iceberg_cmap(round((str2num(datepairs(k).name(1:4))+str2num(datepairs(k).name(16:19)))/2)-2012,:)); hold on;
    plot(Asub,dVdt./86400,'+','markersize',20,'linewidth',2,'markeredgecolor',iceberg_cmap(round((str2num(datepairs(k).name(end-49:end-46))+str2num(datepairs(k).name(end-34:end-31)))/2)-2012,:)); hold on;
%     set(gca,'xlim',[0 max(max_Asub)+0.01*max(max_Asub)]); set(gca,'ylim',[0 max(max_dVdt)+0.01*max(max_dVdt)]); 
    xlims = get(gca,'xlim'); set(gca,'xlim',[0 max(xlims)]);
    xticks = get(gca,'xtick'); set(gca,'xticklabel',xticks/10^6);
    set(gca,'fontsize',20); xlabel('Submerged area (km^2)','fontsize',20); ylabel('Meltwater flux (m^3/s)','fontsize',20);
    grid on;
    subplot(sub2);
%     plot(H,m,'+','markersize',20,'linewidth',2,'markeredgecolor',iceberg_cmap(round((str2num(datepairs(k).name(1:4))+str2num(datepairs(k).name(16:19)))/2)-2012,:)); hold on;
    plot(H,m,'+','markersize',20,'linewidth',2,'markeredgecolor',iceberg_cmap(round((str2num(datepairs(k).name(end-49:end-46))+str2num(datepairs(k).name(end-34:end-31)))/2)-2012,:)); hold on;
%     set(gca,'xlim',[0 max(max_H)+0.01*max(max_H)]); set(gca,'ylim',[0 max(max_melt)+0.01*max(max_melt)]); 
    xlims = get(gca,'xlim'); set(gca,'xlim',[0 max(xlims)]);
    set(gca,'fontsize',20); xlabel('Average iceberg thickness (m)','fontsize',20); ylabel('Melt rate (m/d)','fontsize',20);
    grid on;
    drawnow;
    
    %resave as a tab-delimited text file
%     clear dt xo yo zo po Vo xf yf zf pf Vf coreg_z* dz* dVdt* draft* Asurf* Asub*;
    column_names = {'TimeSeparation (days)' 'X_i (m)' 'Y_i (m)' 'MedianZ_i (m)'...
        'Density_i (kg/m^3)' 'Volume_i (m^3)' 'X_f (m)' 'Y_f (m)' 'MedianZ_f (m)'...
        'Density_f (kg/m^3)' 'Volume_f (m^3)' 'VerticalAdjustment_i (m)' 'VerticalAdjustment_f (m)'...
        'MeanElevationChange (m)' 'StdevElevationChange (m)' 'VolumeChangeRate (m^3/d)' 'VolumeChangeRate_uncert (m^3/d)'...
        'MedianDraft_mean (m)' 'MedianDraft_range (m)' 'SurfaceArea_mean (m^2)' 'SurfaceArea_range (m^2)'...
        'SubmergedArea_mean (m^3)','SubmergedArea_range (m^3)'};
%     append_ref = 1;
%     for i = 1:length(SL)
%         if SL(i).mean.dVdt > 0 && ~isempty(SL(i).mean.TA)
%             dt(append_ref) = sum(SL(i).days);
%             xo(append_ref) = nanmean(SL(i).initial.x); yo(append_ref) = nanmean(SL(i).initial.y); zo(append_ref) = SL(i).initial.z_median; Vo(append_ref) = SL(i).initial.V;
%             xf(append_ref) = nanmean(SL(i).final.x); yf(append_ref) = nanmean(SL(i).final.y); zf(append_ref) = SL(i).final.z_median; Vf(append_ref) = SL(i).final.V;
%             po(append_ref)=SL(i).initial.density; pf(append_ref) = SL(i).final.density;
%             coreg_zo(append_ref) = SL(i).initial.coreg_z; coreg_zf(append_ref) = SL(i).final.coreg_z;
%             dz(append_ref) = SL(i).mean.dz; dz_sigma(append_ref) = SL(i).uncert.dz;
%             dVdt(append_ref) = SL(i).mean.dVdt; dVdt_uncert(append_ref) = max(SL(i).uncert.dVdt);
%             draft(append_ref) = SL(i).mean.draft; draft_uncert(append_ref) = SL(i).change.draft/2;
%             Asurf(append_ref) = SL(i).mean.SA; Asurf_uncert(append_ref) = SL(i).change.SA/2;
%             Asub(append_ref) = SL(i).mean.TA; Asub_uncert(append_ref) = SL(i).change.TA/2;
%             append_ref = append_ref+1;
%         else
%             disp(['Skipping over data for ',num2str(i)]);
%         end
%     end
%     column_vals = [dt' xo' yo' zo' po' Vo' xf' yf' zf' pf' Vf' coreg_zo' coreg_zf' dz' dz_sigma' dVdt' dVdt_uncert' draft' draft_uncert' Asurf' Asurf_uncert' Asub' Asub_uncert'];
    T=table(dt,xo,yo,zo,zo,Vo,xf,yf,zf,pf,Vf,coreg_zo,coreg_zf,dz,dz_sigma,dVdt,dVdt_uncert,draft,2*draft_uncert,Asurf,2*Asurf_uncert,Asub,2*Asub_uncert); T.Properties.VariableNames = column_names;
    
    
    %save in the region directory
%     cd ../..
%     dlmwrite([region_name,'_',datepairs(k).name,'_iceberg_meltinfo.txt'],column_vals,'delimiter','\t');
    writetable(T,[region_name,'_',datepairs(k).name(1:end-3),'csv']);
    
    clear SL dt xo yo zo po Vo xf yf zf pf Vf coreg_zo coreg_zf dz dz_sigma dVdt dVdt_uncert draft draft_uncert Asurf Asurf_uncert Asub Asub_uncert;
end
subplot(sub1); leg = legend(pl,num2str([2013:1:2020]')); set(leg,'location','northwest');
ylims = get(gca,'ylim'); set(gca,'ylim',[0 max(ylims)]);
saveas(gcf,[region_name,'_iceberg-melt_subplots.png'],'png');