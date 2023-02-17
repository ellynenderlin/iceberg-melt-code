function [b,c,vertex_dist,vertex_ang,berg_dist,berg_xoffset,berg_yoffset,k] = measure_iceberg_motion(DEM1_date,DEM2_date,IM1,IM2,A,B,S,cmin,cmax,DEMo_pos,DEMf_pos,imo_pos,imf_pos)
% Function to estimate iceberg motion between DEM times
% Ellyn Enderlin and Rainey Aberle
% Last edit: 17 Feb. 2023
%
% INPUTS:   DEM1_date       DEM1 mean image capture date
%           DEM2_date       DEM2 mean image capture date
%           IM1             structure variable containing earlier
%                               orthoimage info
%           IM2             structure variable containing later orthoimage
%                               info
%           A               DEM1 cropped to iceberg region 
%           B               DEM2 cropped to iceberg region
%           S               iceberg shapefile
%           cmin            minimum elevation in both DEMs
%           cmax            maximum elevation in both DEMs
%           p               iceberg index (for IB structure variable)
%           DEMo_pos        position of figure for plotting DEM1
%           DEMf_pos        position of figure for plotting DEM1
%           imo_pos         position of figure for plotting IM1
%           imf_pos         position of figure for plotting IM2
%
% OUTPUTS:  b               user-determined feature used as the center of 
%                               iceberg rotation in DEM1, visible in both DEMs
%           c               user-determined feature used as the center of 
%                               iceberg rotation in DEM2, visible in both DEMs
%           vertex_dist     distance between each vertex (b and c)
%           vertex_ang      angle between each vertex (b and c)
%           berg_dist       distance between the selected feature (b and c)
%           berg_xoffset    x-direction offset between the selected feature
%           berg_yoffset    y-direction offset between the selected feature
%           k               rotation angle (degrees)
elev_cmap = cmocean('thermal',10001); elev_cmap(1,:) = [1 1 1]; 


%plot DEMs & images to guide polygon rotation
figure1 = figure; set(figure1,'position',DEMo_pos); %set(figure1,'position',[50 800 800 700]);
imagesc(A.x,A.y,A.z_local_adjust); hold on; axis xy equal; set(gca,'clim',[0 cmax],'fontsize',14);
colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
plot(S.X,S.Y,'-*k','linewidth',2,'markersize',4); hold on;
[cont,conth] = contour(A.x,A.y,A.z_local_adjust,[0:1:round(cmax)]);
conth.LineColor = 'k';
title(['Early date: ',num2str(DEM1_date)],'fontsize',16);
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
figureA = figure; set(figureA,'position',imo_pos); %set(figureA,'position',[250 50 400 300]);
imagesc(IM1.x,IM1.y,IM1.z); axis xy equal; colormap gray; hold on;
plot(S.X,S.Y,'-*r','linewidth',2,'markersize',4); hold on;
set(gca,'xlim',xlims,'ylim',ylims);
title(['Early date: ',num2str(DEM1_date)],'fontsize',16);
figure2 = figure; set(figure2,'position',DEMf_pos); %set(figure2,'position',[850 800 800 700]);
imagesc(B.x,B.y,B.z_local_adjust); hold on; axis xy equal; set(gca,'clim',[0 cmax],'fontsize',14);
colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
[cont,conth] = contour(B.x,B.y,B.z_local_adjust,[0:1:round(cmax)]);
conth.LineColor = 'k';
title(['Late date: ',num2str(DEM2_date)],'fontsize',16);
xlims = get(gca,'xlim'); ylims = get(gca,'ylim');
figureB = figure; set(figureB,'position',imf_pos); %set(figureB,'position',[1050 50 400 300]);
imagesc(IM2.x,IM2.y,IM2.z); axis xy equal; colormap gray; hold on;
set(gca,'xlim',xlims,'ylim',ylims);
title(['Late date: ',num2str(DEM2_date)],'fontsize',16);
figure(figure1); figure(figure2);

%pick the axis of rotation
clear b c;
disp('Click on the same (obvious & unique) feature in both DEMs to act as the center for polygon rotation');
disp('...click on the early DEM to activate, then click on the feature');
figure(figure1); ax = waitforbuttonpress; [b] = ginput(1); clear ax;
disp('...click on the later DEM to activate, then click on the feature');
figure(figure2); ax = waitforbuttonpress; [c] = ginput(1); clear ax;

%calculate the distance between each vertex & the selected feature in the earlier DEM
disp('Calculating iceberg translation');
vertex_dist = ((S.X-b(1,1)).^2 + (S.Y-b(1,2)).^2).^(1/2);
vertex_xoffset = S.X - b(1,1); vertex_yoffset = S.Y - b(1,2);
vertex_ang = NaN(1,size(S.X,2));
i = 1;
for i = 1:size(S.X,2)
    if vertex_xoffset(i) > 0 && vertex_yoffset(i) > 0
        vertex_ang(i) = atand(vertex_yoffset(i)/vertex_xoffset(i));
    elseif vertex_xoffset(i) < 0 && vertex_yoffset(i) > 0
        vertex_ang(i) = 180-atand(abs(vertex_yoffset(i)/vertex_xoffset(i)));
    elseif vertex_xoffset(i) < 0 && vertex_yoffset(i) < 0
        vertex_ang(i) = 180+atand(abs(vertex_yoffset(i)/vertex_xoffset(i)));
    else
        vertex_ang(i) = 360-atand(abs(vertex_yoffset(i)/vertex_xoffset(i)));
    end
end

%calculate the distance between the selected feature
berg_dist = ((c(1,1)-b(1,1)).^2 + (c(1,2)-b(1,2)).^2).^(1/2);
berg_xoffset = c(1,1) - b(1,1); berg_yoffset = c(1,2) - b(1,2);

%plot iceberg rotation guiding lines then take an initial guess at how much
%the iceberg rotated between DEM acquisition dates
disp('Calculating iceberg rotation');
% if p==1
    %plot faint berg outlines incrementally rotated on later DEM
    figure(figure2);
    if size(vertex_ang,2) == 1
        vertex_ang = vertex_ang';
    end
    k = [-150:60:150]; %rotate around 0
    for j = 1:length(k)
        new_ang(j,:) = vertex_ang - k(j);
        X_rot(j,:) = vertex_dist.*cosd(new_ang(j,:))+c(1,1);
        Y_rot(j,:) = vertex_dist.*sind(new_ang(j,:))+c(1,2);
        plot(X_rot(j,:),Y_rot(j,:),'linewidth',1,...
            'color',[(abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k))))]); hold on;
        text(X_rot(j,1)+20,Y_rot(j,1),[num2str(k(j)),char(176)],'fontsize',14,...
            'color',[(abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k))))]);
    end
    plot(S.X+berg_xoffset,S.Y+berg_yoffset,'linewidth',2,...
        'color','w'); hold on;
    text(S.X(1)+berg_xoffset+20,S.Y(1)+berg_yoffset,'translated but not rotated',...
        'color','w','fontsize',14);
    clear new_ang X_rot Y_rot;
    %     disp('Take a first guess at how much the iceberg rotated: specify a rotation range by entering ''k = [min:increment:max]; dbcont'' ');
    %     keyboard
    
    
    %determine the best rotation range
    disp('Take a first guess at home much the iceberg rotated between dates (you''ll do this repeatedly!)');
    q=1;
    while q
        rot_prompt = 'Iceberg rotation guess: ';
        k = input(rot_prompt);
        close(figure2);
        
        figure2 = figure;
        imagesc(B.x,B.y,B.z_local_adjust); hold on; axis xy equal; set(gca,'clim',[cmin cmax],'fontsize',14);
        colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
        [cont,conth] = contour(B.x,B.y,B.z_local_adjust,[0:1:round(cmax)]);
        conth.LineColor = 'k';
        title(['Late date: ',num2str(DEM2_date)],'fontsize',16);
        set(figure2,'position',DEMf_pos); %set(figure2,'position',[850 800 800 700]);
        %     set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
        for j = 1:length(k)
            new_ang(j,:) = vertex_ang' - k(j);
            X_rot(j,:) = vertex_dist.*cosd(new_ang(j,:))+c(1,1);
            Y_rot(j,:) = vertex_dist.*sind(new_ang(j,:))+c(1,2);
            plot(X_rot(j,:),Y_rot(j,:),'linewidth',2,...
                'color',[0.5 0.5 0.5]); hold on;
            if j == 1 || j == length(k)
                text(X_rot(j,1)+20,Y_rot(j,1),[num2str(k(j)),char(176)],'fontsize',14,...
                    'color',[0.5 0.5 0.5]);
            end
        end
        clear new_ang X_rot Y_rot;
        
        prompt = 'Are you happy with the rotation estimate (y/n)?'; str = input(prompt,'s');
        if strmatch(str,'n')==1
            q = q+1;
        else
            break
        end
    end
    
    
    
    %     rot_prompt = 'Take a first guess at how much the iceberg rotated ([min:increment:max]): ';
    %     k = input(rot_prompt);
    %     close(figure2);
    %
    %     %specify a narrower range of rotation by comparing rotated polygons
    %     figure2 = figure;
    %     imagesc(B.x,B.y,B.z_local_adjust); hold on; axis xy equal; set(gca,'clim',[cmin cmax],'fontsize',14);
    %     colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
    %     [cont,conth] = contour(B.x,B.y,B.z_local_adjust,[0:1:round(cmax)]);
    %     conth.LineColor = 'k';
    %     title(['Late date: ',num2str(DEM2_date)],'fontsize',16);
    %     set(figure2,'position',DEMf_pos); %set(figure2,'position',[850 800 800 700]);
    %     %     set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
    %     for j = 1:length(k)
    %         new_ang(j,:) = vertex_ang' - k(j);
    %         X_rot(j,:) = vertex_dist.*cosd(new_ang(j,:))+c(1,1);
    %         Y_rot(j,:) = vertex_dist.*sind(new_ang(j,:))+c(1,2);
    %         plot(X_rot(j,:),Y_rot(j,:),'linewidth',2,...
    %             'color',[(abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k))))]); hold on;
    %         if j == 1 || j == length(k)
    %             text(X_rot(j,1)+20,Y_rot(j,1),[num2str(k(j)),char(176)],'fontsize',14,...
    %                 'color',[(abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k))))]);
    %         end
    %     end
    %     clear new_ang X_rot Y_rot;
    
% end

% %narrow the range of rotation again by comparing rotated polygons
% % disp('Set a narrower rotation range by entering ''k = [min:increment:max]; dbcont'' ');
% % keyboard
% rot_prompt = 'Set a narrower rotation range ([min:increment:max]): ';
% k = input(rot_prompt);
% close(figure2);
% figure2 = figure; set(figure2,'position',DEMf_pos); %set(figure2,'position',[850 800 800 700]);
% imagesc(B.x,B.y,B.z_local_adjust); hold on; axis xy equal; set(gca,'clim',[cmin cmax],'fontsize',14);
% colormap(gca,elev_cmap); cbar = colorbar; set(get(cbar,'ylabel'),'string', 'elevation (m)');
% [cont,conth] = contour(B.x,B.y,B.z_local_adjust,[0:1:round(cmax)]);
% conth.LineColor = 'k';
% title(['Late date: ',num2str(DEM2_date)],'fontsize',16);
% % set(gca,'xlim',[min(a(:,1)) max(a(:,1))],'ylim',[min(a(:,2)) max(a(:,2))]);
% for j = 1:length(k)
%     new_ang(j,:) = vertex_ang' - k(j);
%     X_rot(j,:) = vertex_dist.*cosd(new_ang(j,:))+c(1,1);
%     Y_rot(j,:) = vertex_dist.*sind(new_ang(j,:))+c(1,2);
%     figure(figure2);
%     plot(X_rot(j,:),Y_rot(j,:),'linewidth',2,...
%         'color',[(abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k))))]); hold on;
%     if j == 1 || j == length(k)
%         text(X_rot(j,1)+20,Y_rot(j,1),[num2str(k(j)),char(176)],'fontsize',14,...
%             'color',[(abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k))))]);
%     end
%     figure(figureB);
%     plot(X_rot(j,:),Y_rot(j,:),'linewidth',2,...
%         'color',[1 (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k))))]); hold on;
%     if j == 1 || j == length(k)
%         text(X_rot(j,1)+20,Y_rot(j,1),[num2str(k(j)),char(176)],'fontsize',14,...
%             'color',[1 (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k)))) (abs(k(j))-(min(abs(k))))/((max(abs(k))-min(abs(k))) + round(0.25*max(abs(k))))]);
%     end
% end
% clear new_ang X_rot Y_rot;
% 
% %pick the final rotation
% % disp('Set the final iceberg rotation estimate by entering ''k = XX; dbcont'' ');
% % keyboard
% rot_prompt = 'Set the final iceberg rotation estimate: ';
% k = input(rot_prompt);
% if length(k) == 1
%     return
% end

end