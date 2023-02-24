function [PSx_early,PSy_early,PSx_late,PSy_late] = track_icebergs(DEM1,DEM2,IM1,IM2,dir_output)
% Function to identify the identifiable icebergs in both DEMs and save the 
% coordinates as .txt files.
% Ellyn Enderlin (ellynenderlin@boisestate.edu)
% Reformatted by Rainey Aberle, Fall 2021
%
% INPUTS:   DEM1            structure variable containing earlier DEM info
%           DEM2            structure variable containing later DEM info
%           IM1             structure variable containing earlier
%                               orthoimage info
%           IM2             structure variable containing later orthoimage
%                               info
%           dir_output      directory where all output files will be placed
%
% OUTPUTS:  PSx_early       x coordinates (Antarctic Polar Stereo) of
%                           icebergs in the earlier DEM
%           PSy_early       y coordinates (Antarctic Polar Stereo) of
%                           icebergs in the earlier DEM
%           PSx_late        x coordinates (Antarctic Polar Stereo) of
%                           icebergs in the later DEM
%           PSy_late        y coordinates (Antarctic Polar Stereo) of
%                           icebergs in the later DEM

%specify reasonable iceberg height above ocean surface for color-scaling
maxelev = 20;
elev_cmap = cmocean('thermal',1001);

%navigate to the date-specific directory
cd_to_datefolder = ['cd ',dir_output,'/',DEM1.time,'-',DEM2.time,'/']; eval(cd_to_datefolder);

%plot co-registered DEMs
% DEM1.z_elpsd_adjust(DEM1.z_elpsd_adjust<0) = 0; DEM2.z_elpsd_adjust(DEM2.z_elpsd_adjust<0) = 0; %remove elevations less than zero for plotting purposes
figure1 = figure; set(figure1,'position',[0 600 800 600]);
imagesc(DEM1.x,DEM1.y,DEM1.z_elpsd_adjust); axis xy equal; grid on; hold on;
set(gca,'clim',[min(nanmedian(DEM1.z_elpsd_adjust)) min(nanmedian(DEM1.z_elpsd_adjust))+maxelev]); colormap(gca,elev_cmap); colorbar; 
set(gca,'xtick',[min(DEM1.x):round(range(DEM1.x)/10000)*1000:max(DEM1.x)],'xticklabel',[min(DEM1.x)/1000:round(range(DEM1.x)/10000):max(DEM1.x)/1000],...
    'ytick',[min(DEM1.y):round(range(DEM1.y)/10000)*1000:max(DEM1.y)],'yticklabel',[min(DEM1.y)/1000:round(range(DEM1.y)/10000):max(DEM1.y)/1000]);
figure2 = figure; set(figure2,'position',[975 600 800 600]);
imagesc(DEM2.x,DEM2.y,DEM2.z_elpsd_adjust); axis xy equal; grid on; hold on;
set(gca,'clim',[min(nanmedian(DEM2.z_elpsd_adjust)) min(nanmedian(DEM2.z_elpsd_adjust))+maxelev]); colormap(gca,elev_cmap); colorbar; 
set(gca,'xtick',[min(DEM2.x):round(range(DEM2.x)/10000)*1000:max(DEM2.x)],'xticklabel',[min(DEM2.x)/1000:round(range(DEM2.x)/10000):max(DEM2.x)/1000],...
    'ytick',[min(DEM2.y):round(range(DEM2.y)/10000)*1000:max(DEM2.y)],'yticklabel',[min(DEM2.y)/1000:round(range(DEM2.y)/10000):max(DEM2.y)/1000]);

%plot images
figure3 = figure; set(figure3,'position',[0 0 800 600]); ax3=gca;
imagesc(IM1.x,IM1.y,IM1.z); axis xy equal; colormap(gca,'gray'); grid on; set(gca,'clim',[0 200]); hold on;
set(gca,'xtick',[min(DEM1.x):round(range(DEM1.x)/10000)*1000:max(DEM1.x)],'xticklabel',[min(DEM1.x)/1000:round(range(DEM1.x)/10000):max(DEM1.x)/1000],...
    'ytick',[min(DEM1.y):round(range(DEM1.y)/10000)*1000:max(DEM1.y)],'yticklabel',[min(DEM1.y)/1000:round(range(DEM1.y)/10000):max(DEM1.y)/1000]);

figure4 = figure; set(figure4,'position',[975 0 800 600]); ax4=gca;
imagesc(IM2.x,IM2.y,IM2.z); axis xy equal; colormap(gca,'gray'); grid on; set(gca,'clim',[0 200]); hold on;
set(gca,'xtick',[min(DEM2.x):round(range(DEM2.x)/10000)*1000:max(DEM2.x)],'xticklabel',[min(DEM2.x)/1000:round(range(DEM2.x)/10000):max(DEM2.x)/1000],...
    'ytick',[min(DEM2.y):round(range(DEM2.y)/10000)*1000:max(DEM2.y)],'yticklabel',[min(DEM2.y)/1000:round(range(DEM2.y)/10000):max(DEM2.y)/1000]);


% pick iceberg coordinates
disp('Follow prompts to Zoom in & out on the images & select icebergs for analysis.')
disp('Try to focus on big icebergs. Max icebergs = 30!');

%check for existing coordinates
PScoord_files = dir([dir_output,DEM1.time,'-',DEM2.time,'/','*PScoords.txt']);
if ~isempty(PScoord_files)
    disp(['Already matched ',num2str(length(PScoord_files)),' icebergs']);
    for q = 1:length(PScoord_files)
        %read in the coordinates
        iceberg_coords = [PScoord_files(q).folder,'/',PScoord_files(q).name];
        coords = cell2mat(textscan(fopen(iceberg_coords),'%f64 %f64 %f64 %f64','Delimiter',',','headerlines',1));
        PSy_early(q) = coords(1); PSx_early(q) = coords(2);
        PSy_late(q) = coords(3); PSx_late(q) = coords(4);

        %plot coordinates
        figure(figure3);
        plot(PSx_early(q),PSy_early(q),'.','markersize',10,'color','c');
        figure(figure1);
        plot(PSx_early(q),PSy_early(q),'.','markersize',10,'color','c');
        drawnow;
        disp('later image');
        figure(figure4);
        plot(PSx_late(q),PSy_late(q),'.','markersize',10,'color','c');
        figure(figure2);
        plot(PSx_late(q),PSy_late(q),'.','markersize',10,'color','c');
        clear coords;
    end
    q = length(PScoord_files)+1;
else
    q=1;
end

%iterative iceberg selection
while q
    ID_question = questdlg('Are there matching icebergs?',...
        'Iceberg Pairing','1) Yes!','2) No!','1) Yes!');
    if q <= 30
        %execute blunder removal based on question response
        switch ID_question
            case '1) Yes!'
                disp(['iceberg #',num2str(q)]);
                
                %zoom in
                disp('earlier image');
                figure(figure3);
                disp('Click on UL & LR corners of a box bounding a region where you want to look at icebergs to zoom in'); % Upper left, lower right.
                [a] = ginput(2);
                set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
                figure(figure1);
                set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
                drawnow; clear a;
                disp('later image');
                figure(figure4);
                disp('Click on UL & LR corners of a box bounding a region where you want to look at icebergs to zoom in'); % Upper left, lower right.
                [a] = ginput(2);
                set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
                figure(figure2);
                set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
                drawnow; clear a;
                
                disp('Click on the middle of the iceberg in each DEM');
                disp('earlier image');
                figure(figure3);
                [PSx_early(q),PSy_early(q)] = ginput(1);
                plot(PSx_early(q),PSy_early(q),'.','markersize',10,'color','c');
                figure(figure1);
                plot(PSx_early(q),PSy_early(q),'.','markersize',10,'color','c');
                drawnow;
                disp('later image');
                figure(figure4);
                [PSx_late(q),PSy_late(q)] = ginput(1);
                plot(PSx_late(q),PSy_late(q),'.','markersize',10,'color','c');
                figure(figure2);
                plot(PSx_late(q),PSy_late(q),'.','markersize',10,'color','c');
                drawnow;
                
                %zoom back out & advance
                q = q+1;
                figure(figure1); set(gca,'xlim',[min(DEM1.x) max(DEM1.x)],'ylim',[min(DEM1.y) max(DEM1.y)]);
                figure(figure3); set(gca,'xlim',[min(DEM1.x) max(DEM1.x)],'ylim',[min(DEM1.y) max(DEM1.y)]);
                figure(figure2); set(gca,'xlim',[min(DEM2.x) max(DEM2.x)],'ylim',[min(DEM2.y) max(DEM2.y)]);
                figure(figure4); set(gca,'xlim',[min(DEM2.x) max(DEM2.x)],'ylim',[min(DEM2.y) max(DEM2.y)]);
                drawnow;
            case '2) No!'
                disp('Move on!')
                break
        end
    else
        disp('Matched 30 icebergs. Move on!');
        break
    end
    clear ID_question;
end

% % col = spring(30); % color scheme for plotting
% more=1; % select more icebergs when more=1
% i=0; % index for storing coordinates
% axes(ax3); % bring figure 3 axes to front
% while more==1
%     axes(ax3);
%     i=i+1;
%     [PSx_early(i),PSy_early(i)] = ginput(1);
%     plot(gca,PSx_early(i),PSy_early(i),'.','markersize',10,'color','y');
%     
%     str = input('Do you see more icebergs (y/n)?','s');
%     if strcmp(str,'n')
%         more=more+1;
%     end 
% end
% %add number labels
% axes(ax3);
% for i = 1:length(PSx_early)
%    text(PSx_early(i)+300,PSy_early(i),num2str(i),'color','y','fontsize',18); hold on;
% end
% disp('later DEM: Now, click on the same icebergs in the other image');
% disp('   You can click until the same number of icebergs have been selected');
% figure(figure4);
% for j=1:i
%     axes(ax4); % bring figure 4 axes to front
%     [PSx_late(j),PSy_late(j)] = ginput(1);
%     plot(PSx_late(j),PSy_late(j),'.','markersize',10,'color','y');
%     drawnow
% end

% concatenate coordinates
early_coords = [PSx_early',PSy_early']; %early_coords = sortrows(early_coords,[1 2]);
PSx_early = early_coords(:,1); PSy_early = early_coords(:,2);
late_coords = [PSx_late',PSy_late']; %late_coords = sortrows(late_coords,[1 2]);
PSx_late = late_coords(:,1); PSy_late = late_coords(:,2);

%save iceberg coordinates
% if ~exist([dir_output,'iceberg_data/'], 'dir')
%    mkdir([dir_output,'iceberg_data/'])
% end
for j = 1:length(PSx_early)
    if size(num2str(j),2) == 1
        iceberg_no = [num2str(0),num2str(j)];
    else
        iceberg_no = num2str(j);
    end
    coords = [PSy_early(j)  PSx_early(j) PSy_late(j) PSx_late(j)];
    coords_table = array2table(coords,...
        'VariableNames',{'DEM1: Y (m)','DEM1: X (m)','DEM2: Y (m)','DEM2: X (m)'});
    writetable(coords_table,[dir_output,'/',DEM1.time,'-',DEM2.time,'/','iceberg',iceberg_no,'_PScoords.txt']);
    clear coords;
end
disp('iceberg coordinates saved');

% close all; drawnow;
clear Y Z I1 I2;

disp('Advance to the next step');
disp('------------------------');