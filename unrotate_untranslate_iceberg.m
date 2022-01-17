function [x_rot,y_rot] = unrotate_untranslate_iceberg(IB,p)
% Function to unrotate and untranslate iceberg between DEM times
% Ellyn Enderlin, Mariama Dryak
% Slightly modified by Rainey Aberle, Fall 2021
%
% INPUTS:   IB      structure variable with iceberg info
%           p       number of image pixels in iceberg polygon
% 
% OUTPUTS:  x_rot   new rotated iceberg x coordinates
%           y_rot   new rotated iceberg y coordinates

%calculate the translation & rotation needed to adjust each pixel to the
%same reference frame
for i = 1:size(IB(p).zf.local_adjust.map,1)
    for k = 1:size(IB(p).zf.local_adjust.map,2)
        dist(i,k) = ((IB(p).xf(1,k)-IB(p).vertices.xf_mid).^2 + (IB(p).yf(1,i)-IB(p).vertices.yf_mid).^2).^(1/2);
        xoffset(1,k) = IB(p).xf(1,k) - IB(p).vertices.xf_mid; yoffset(1,i) = IB(p).yf(1,i) - IB(p).vertices.yf_mid;
        ang(i,k) = NaN;
    end
end
%calculate the angle
for i = 1:size(IB(p).zf.local_adjust.map,1);
    for k = 1:size(IB(p).zf.local_adjust.map,2)
        if xoffset(k) > 0 && yoffset(i) > 0;
            ang(i,k) = atand(yoffset(i)/xoffset(k));
        elseif xoffset(k) < 0 && yoffset(i) > 0;
            ang(i,k) = 180-atand(abs(yoffset(i)/xoffset(k)));
        elseif xoffset(k) < 0 && yoffset(i) < 0;
            ang(i,k) = 180+atand(abs(yoffset(i)/xoffset(k)));
        else
            ang(i,k) = 360-atand(abs(yoffset(i)/xoffset(k)));
        end
    end
end

%rotate and shift the iceberg from the later DEM to the earlier DEM
%coordinates
new_ang = ang - IB(p).rotate; %correct for rotation
for f = 1:size(IB(p).zf.local_adjust.map,1);
    for g = 1:size(IB(p).zf.local_adjust.map,2)
        if new_ang(f,g) >= 0 && new_ang(f,g) < 90;
            x_rot(f,g) = dist(f,g)*cosd(new_ang(f,g)) + IB(p).vertices.xo_mid;
            y_rot(f,g) = dist(f,g)*sind(new_ang(f,g)) + IB(p).vertices.yo_mid;
        elseif new_ang(f,g) >= 90 && new_ang(f,g) < 180;
            x_rot(f,g) = -dist(f,g)*cosd(180-new_ang(f,g)) + IB(p).vertices.xo_mid;
            y_rot(f,g) = dist(f,g)*sind(180-new_ang(f,g)) + IB(p).vertices.yo_mid;
        elseif new_ang(f,g) >= 180 && new_ang(f,g) < 270;
            x_rot(f,g) = -dist(f,g)*cosd(new_ang(f,g)-180) + IB(p).vertices.xo_mid;
            y_rot(f,g) = -dist(f,g)*sind(new_ang(f,g)-180) + IB(p).vertices.yo_mid;
        else
            x_rot(f,g) = dist(f,g)*cosd(360-new_ang(f,g)) + IB(p).vertices.xo_mid;
            y_rot(f,g) = -dist(f,g)*sind(360-new_ang(f,g)) + IB(p).vertices.yo_mid;
        end
    end
end

end