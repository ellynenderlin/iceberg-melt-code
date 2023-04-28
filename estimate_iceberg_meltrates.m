function [DEM1,DEM2] = estimate_iceberg_meltrates(DEM1,DEM2,IM1,IM2,dir_output,dir_code,dir_SMB,geography,region_abbrev,region_name,option_no)
% Function to estimate iceberg freshwater fluxes and melt rates
% Ellyn Enderlin (ellynenderlin@boisestate.edu) and Rainey Aberle (raineyaberle@u.boisestate.edu)
% Last edited: 09 Nov. 2022
%
% INPUTS:   DEM1            structure variable containing earlier DEM info
%           DEM2            structure variable containing later DEM info
%           IM1             structure variable containing earlier
%                               orthoimage info
%           IM2             structure variable containing later orthoimage
%                               info
%           dir_output      directory where all output files will be placed
%           dir_code        directory to the Iceberg-melt-rate code folder
%                               (including the name of the folder)
%           geography       binary specification of polar region (0 = Greenland, 1 = Antarctic)
%           region_abbrev   region abbreviation used in image files
%           option_no         select which step to execute (1 or 2):
%                           (1) Estimate elevation change for each iceberg
%                           (2) Update individual icebergs and/or remove
%                               icebergs with anomalous melt rate estimates
%
% OUTPUTS:  DEM1            structure variable containing earlier DEM info
%                               with any new fields
%           DEM2            structure variable containing later DEM info
%                               with any new fields
%
% Calls the following external functions:
%   - extract_iceberg_elev_change.m
%   - convert_iceberg_elev_change_to_meltrates.m


% ----------STEP 1: Calculate Elevation Change----------
dir_iceberg = [dir_output,DEM1.time,'-',DEM2.time,'/'];
cd(dir_iceberg);

%select the iceberg number for which to start the elevation change estimates
if option_no ~= 3
    disp('Extract iceberg elevation change');
    icebergs = dir([dir_iceberg,'iceberg*coords.txt']);
    iceberg_dz = dir([dir_iceberg,'iceberg*dz.mat']);
    if ~isempty(iceberg_dz)
        disp(['Already calculated elevation change for ',num2str(length(iceberg_dz)),' of ',num2str(length(icebergs)),' icebergs']);
        answer = questdlg('Do you want/need to calculate elevation changes for more icebergs?',...
            'Elevation Change Estimation','1) Yes!','2) No','1) Yes!');
        switch answer
            case '1) Yes!'
                %specify the iceberg numbers to loop through
                disp('Specify range of iceberg numbers as "iceberg_refs = X:Y; dbcont" in the command window (w/o quotes) then hit enter to loop');
                disp('   Ex of loop (for 4 through the rest): iceberg_refs = 4:size(icebergs,1); dbcont');
                disp('   Ex of select numbers: iceberg_refs = [5,15,16]; dbcont');
                keyboard

                %loop
                for j = iceberg_refs %size(icebergs,1) %default to loop through all icebergs is j = 1:size(icebergs,1)
                    if j<10; iceberg_no = ['0',num2str(j)]; else iceberg_no = num2str(j); end
                    [IB,dz] = extract_iceberg_elev_change(DEM1,DEM2,IM1,IM2,iceberg_no,dir_output,dir_code,geography,region_abbrev);
                    clear IB dz;
                end
            case '2) No'
                disp('Moving on to melt rate estimation...');
        end
    else
        for j = 1:size(icebergs,1) %size(icebergs,1) %default to loop through all icebergs is j = 1:size(icebergs,1)
            iceberg_no = icebergs(j).name(8:9);
            [IB,dz] = extract_iceberg_elev_change(DEM1,DEM2,IM1,IM2,iceberg_no,dir_output,dir_code,geography,region_abbrev);
            clear IB dz;
        end
    end
    close all;
end

% ----------STEP 2: Estimate Melt Rates----------
if option_no==1 %calculate melt rates for all icebergs
    
    % Calculate meltwater fluxes, melt rates, and uncertainties
    disp('Convert elevation change to meltwater fluxes & melt rates');
    cd(dir_iceberg);
    berg_numbers = dir([dir_iceberg,'iceberg*dz.mat']);
    dir_bedrock = [dir_output,'DEM_offset_data/'];
    
    % calculate iceberg meltrate from elevation change:
    if exist([region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat']) == 2 %does the meltrate file already exist?
        redo = questdlg('Do you want/need to redo the conversion to melt rates?',...
            'Redo melt conversion','1) Yes','2) No','1) Yes');
        switch redo
            case '1) Yes'
                SL = convert_iceberg_elev_change_to_meltrates(DEM1,DEM2,IM1,IM2,berg_numbers,geography,region_name,region_abbrev,dir_output,dir_code,dir_iceberg,dir_SMB);
            case '2) No'
                disp('reloading meltrate data...');
                load([region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat']);
        end
    else
        SL = convert_iceberg_elev_change_to_meltrates(DEM1,DEM2,IM1,IM2,berg_numbers,geography,region_name,region_abbrev,dir_output,dir_code,dir_iceberg,dir_SMB);
    end
    
    %plot & export to table
    plot_flag = 1; %plot data
    table_flag = 1; %export summary data to table
    plot_export_iceberg_melt_data(SL,dir_output,region_abbrev,DEM1,DEM2,plot_flag,table_flag);
    
    %call-out the clearly bad icebergs
    for i = 1:length(SL)
        if SL(i).mean.dVdt < 0
            disp(['Recalculate elevation change for iceberg',SL(i).name(end-1:end)]);
        end
    end
    
    %resave to the mat-file
    disp('Saving melt rates');
    save([dir_output,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat'],'SL','-v7.3');
    
    disp(' ');
    disp('Run option 2 in the last section of the wrapper:');
    disp('   a) When prompted, specify the icebergs to rerun using iceberg_refs = X:Y; dbcont');
    disp('   b) Update and/or remove icebergs as necessary');
    clear SL;
    
elseif option_no==2 %recalculate melt rates for select icebergs
    disp('Recalculate melt rates for select icebergs then remove icebergs that still have bad results');
    dir_iceberg = [dir_output,DEM1.time,'-',DEM2.time,'/'];
    
    %update and/or remove select icebergs
    if exist('iceberg_refs') ~= 1
        disp('Specify the berg numbers (from melt plots) that need updating as "iceberg_refs = []; dbcont"');
        keyboard
    end
    SL = update_iceberg_meltrates(DEM1,DEM2,IM1,IM2,geography,region_name,region_abbrev,iceberg_refs,dir_output,dir_iceberg,dir_code,dir_SMB);
    disp('Resaved iceberg melt structure, plots, & data table!');

elseif option_no==3 %only plot
    load([dir_output,region_abbrev,'_',DEM1.time,'-',DEM2.time,'_iceberg_melt.mat']);
    cd(dir_iceberg);
    plot_flag = 1; %plot data
    table_flag = 0; %default (0) suppresses table generation, switch to 1 to create tables
    plot_export_iceberg_melt_data(SL,dir_output,region_abbrev,DEM1,DEM2,plot_flag,table_flag);
end

end