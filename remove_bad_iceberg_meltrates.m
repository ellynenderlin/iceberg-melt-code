function [SL] = remove_bad_iceberg_meltrates(DEM1,DEM2,region_abbrev,iceberg_refs,dir_output)
% Function to remove icebergs with melt rates clearly biased by poor sea
% level adjustments based on inspection of melt rate vs. draft and melt
% rate vs. sea level correction plots.
% Ellyn Enderlin, Summer 2022
%
% INPUTS:   DEM1            structure variable containing earlier DEM info
%           DEM2            structure variable containing later DEM info
%           iceberg_refs    number of iceberg to detect elevation change
%           dir_output      directory where all output files will be placed
%           region_abbrev   abbrevation of region used in file names
%
% OUTPUTS:  Updates the SL structure containing iceberg melt volume & melt
% rate data to remove flagged icebergs

%load the iceberg melt rate (SL) structure
DEM1_date = DEM1.time(1:8); DEM2_date = DEM2.time(1:8);
load([dir_output,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat']);

%identify the structure reference for the specified icebergs
for i = 1:length(iceberg_refs)
    for k = 1:length(SL)
        if iceberg_refs(i) < 10
            if strmatch([num2str(0),num2str(iceberg_refs(i))],SL(k).name(end-1:end))
                berg_refs(i) = k;
            end
        else
            if strmatch(num2str(iceberg_refs(i)),SL(k).name(end-1:end))
                berg_refs(i) = k;
            end
        end
    end
end

%remove the bad icebergs from the structure
for i = 1:length(berg_refs)
    SL(berg_refs(i)).mean.TA = [];
    SL(berg_refs(i)).mean.dHdt = [];
end

%resave without the bad icebergs
save([dir_output,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat'],'SL','-v7.3');

%replot and re-export data to a table
table_flag = 1; %save summary data to table
plot_flag = 1; %plot figures
[T] = plot_export_iceberg_melt_data(SL,dir_output,region_abbrev,DEM1,DEM2,plot_flag,table_flag);
fprintf(['Iceberg melt data saved to table from ',region_abbrev,' for %i icebergs \n'],height(T));


end

