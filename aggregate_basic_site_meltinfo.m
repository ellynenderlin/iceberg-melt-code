%append more data to iceberg melt csvs
clearvars; close all;

%specify directories for input CSVs produced by the iceberg melt mapping
%pipeline & for the site-concatenated CSVs
input_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/iceberg-melt-data/original-csvs/';
output_dir = '/Users/ellynenderlin/Research/NSF_GrIS-Freshwater/iceberg-melt-data/site-csvs/';

%identify the meltrate csv files & extract site names
csvs = dir([input_dir,'*iceberg_meltinfo.csv']);
for j = 1:length(csvs)
    site_abbrev(j,:) = csvs(j).name(1:3);
end

%find file indices for each unique site name
[unique_abbrevs,ia,ix] = unique(site_abbrev,'rows');

%create site-specific CSVs with additional columns of data
for j = 1:length(ia)
    disp(unique_abbrevs(j,:));
    site_refs = find(ix == j);
    T = table;
    
    %loop through all the dated CSVs for the site
    for k = 1:length(site_refs)
        %read the table
        disp(csvs(site_refs(k)).name);
        Ttemp = readtable([input_dir,csvs(site_refs(k)).name]);
        
        %make a new table with less info
        Td = table;
        
        Td.('Site Abbreviation') = repmat(unique_abbrevs(j,:),size(Ttemp.VolumeChangeRate));
        Td.('Date start (YYYYMMDD)') = repmat(csvs(site_refs(k)).name(5:12),size(Ttemp.VolumeChangeRate));
        Td.('Date end (YYYYMMDD)') = repmat(csvs(site_refs(k)).name(14:21),size(Ttemp.VolumeChangeRate));
        Td.('X (m)') = mean([Ttemp.X_i,Ttemp.X_f],2);
        Td.('Y (m)') = mean([Ttemp.Y_i,Ttemp.Y_f],2);
        Td.('Z (m)') = mean([Ttemp.MedianZ_i,Ttemp.MedianZ_f],2);
        Td.('Density (kg m^-3)') = mean([Ttemp.Density_i,Ttemp.Density_f],2);
        Td.('Volume start (m^3)') = Ttemp.Volume_i;
        Td.('Volume end (m^3)') = Ttemp.Volume_f;
        Td.('Vertical coreg. start (m)') = Ttemp.VerticalAdjustment_i;
        Td.('Vertical coreg. end (m)') = Ttemp.VerticalAdjustment_f;
        Td.('dZ mean (m)') = Ttemp.ElevationChange_mean;
        Td.('dZ st.dev. (m)') = Ttemp.ElevationChange_stdev;
        Td.('dVdt mean (m^3 d^-1)') = Ttemp.VolumeChangeRate;
        Td.('dVdt uncert. (m^3 d^-1)') = Ttemp.VolumeChangeRate_uncert;
        Td.('Draft mean (m)') = Ttemp.MedianDraft_mean;
        Td.('Draft range (m)') = Ttemp.MedianDraft_range;
        Td.('Surface Area mean (m^2)') = Ttemp.SurfaceArea_mean;
        Td.('Surface Area range (m^2)') = Ttemp.SurfaceArea_range;
        Td.('Submerged Area mean (m^2)') = Ttemp.SubmergedArea_mean;
        Td.('Submerged Area uncert. (m^2)') = Ttemp.SubmergedArea_uncert;
        
        %calculate melt rate (m/d)
        Td.('MeltRate (m d^-1)') = Ttemp.VolumeChangeRate./Ttemp.SubmergedArea_mean;
        Td.('MeltRate uncert. (m d^-1)') = abs(Ttemp.VolumeChangeRate./Ttemp.SubmergedArea_mean).*sqrt((Ttemp.VolumeChangeRate_uncert./Ttemp.VolumeChangeRate).^2 + (Ttemp.SubmergedArea_uncert./Ttemp.SubmergedArea_mean).^2);
        
        %add to the site table
        T = [T;Td];
        
        clear Ttemp Td;
    end
    
    %export the new table as a CSV
    writetable(T,[output_dir,unique_abbrevs(j,:),'_basic_iceberg_meltinfo.csv']);
    clear T;
end
disp('Done creating site-specific iceberg melt files!');


