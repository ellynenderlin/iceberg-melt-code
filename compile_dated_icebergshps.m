%%% Write a single shapefile for each site datepair
clearvars; close all;

%testing for a single datepair at a single site before iterating
cd /Users/Shared/Greenland/melange/

%iterate through sites, looking for datepair subdirectories
sites = dir('*');
for p = 1:length(sites)
    if length(sites(p).name) == 3 && isfolder(sites(p).name)
        cd(sites(p).name); disp(sites(p).name);

        %iterate through datepair subdirectories
        datepairs = dir('20*-20*');
        for k = 1:length(datepairs)
            cd(datepairs(k).name);

            %check for the shapefiles directory w/ shapes it in
            if isfolder('iceberg_shapes')
                cd iceberg_shapes

                %iterate through the shapefiles & combine
                shps = dir('WV*.shp'); S = [];
                if ~isempty(shps)
                    for j = 1:length(shps)
                        Stemp = readgeotable(shps(j).name);
                        Stemp.SiteID = sites(p).name; %add the 3-letter site id
                        Stemp.Date = shps(j).name(4:11); %add the observation date
                        Stemp.IcebergID = shps(j).name(end-5:end-4); %add the iceberg id
                        Stemp.Name = []; %remove the Name field ("iceberg##")
                        S = [S; Stemp]; % Append the temporary table to the main table

                        clear Stemp;
                    end
                    shapewrite(S,[sites(p).name,'_',shps(1).name(4:11),'-',shps(end).name(4:11),'_icebergs.shp']);
                    copyfile([shps(1).name(1:end-3),'prj'],[sites(p).name,'_',shps(1).name(4:11),'-',shps(end).name(4:11),'_icebergs','.prj']);
                    disp(['file written for ',shps(1).name(4:11),'-',shps(end).name(4:11)]);
                end
                clear S shps;

                cd ..
            end

            cd ..
        end
        clear datepairs;

        cd ..
    end

end

