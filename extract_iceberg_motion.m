function [b,c,vertex_dist,vertex_ang,berg_dist,berg_xoffset,berg_yoffset,k] = extract_iceberg_motion(S,Sl)

%calculate the distance between each vertex & the selected feature in the
%earlier DEM
disp('Extracting iceberg polygon vertex information');
vertex_dist = ((S.X-nanmean(S.X)).^2 + (S.Y-nanmean(S.Y)).^2).^(1/2);
vertex_xoffset = S.X - nanmean(S.X); vertex_yoffset = S.Y - nanmean(S.Y);
vertex_ang = NaN(1,size(S.X,2));
i = 1;
for i = 1:size(S.X,2);
    if vertex_xoffset(i) > 0 && vertex_yoffset(i) > 0;
        vertex_ang(i) = atand(vertex_yoffset(i)/vertex_xoffset(i));
    elseif vertex_xoffset(i) < 0 && vertex_yoffset(i) > 0;
        vertex_ang(i) = 180-atand(abs(vertex_yoffset(i)/vertex_xoffset(i)));
    elseif vertex_xoffset(i) < 0 && vertex_yoffset(i) < 0;
        vertex_ang(i) = 180+atand(abs(vertex_yoffset(i)/vertex_xoffset(i)));
    else
        vertex_ang(i) = 360-atand(abs(vertex_yoffset(i)/vertex_xoffset(i)));
    end
end

%read later DEM shape-file iceberg polygon vertex info 
vertex_xoffset_l = Sl.X - nanmean(Sl.X); vertex_yoffset_l = Sl.Y - nanmean(Sl.Y);
vertex_ang_l = NaN(1,size(Sl.X,2));
i = 1;
for i = 1:size(Sl.X,2);
    if vertex_xoffset_l(i) > 0 && vertex_yoffset_l(i) > 0;
        vertex_ang_l(i) = atand(vertex_yoffset_l(i)/vertex_xoffset_l(i));
    elseif vertex_xoffset_l(i) < 0 && vertex_yoffset_l(i) > 0;
        vertex_ang_l(i) = 180-atand(abs(vertex_yoffset_l(i)/vertex_xoffset_l(i)));
    elseif vertex_xoffset_l(i) < 0 && vertex_yoffset_l(i) < 0;
        vertex_ang_l(i) = 180+atand(abs(vertex_yoffset_l(i)/vertex_xoffset_l(i)));
    else
        vertex_ang_l(i) = 360-atand(abs(vertex_yoffset_l(i)/vertex_xoffset_l(i)));
    end
end

%calculate the distance between the selected feature
berg_dist = ((nanmean(Sl.X)-nanmean(S.X)).^2 + (nanmean(Sl.Y)-nanmean(S.Y)).^2).^(1/2);
berg_xoffset = nanmean(Sl.X) - nanmean(S.X); berg_yoffset = nanmean(Sl.Y) - nanmean(S.Y);
rotate_ang = vertex_ang - vertex_ang_l;
k=roundn(nanmean(rotate_ang),-1);

b(1,1) = nanmean(S.X); b(1,2) = nanmean(S.Y);
c(1,1) = nanmean(Sl.X); c(1,2) = nanmean(Sl.Y);

end