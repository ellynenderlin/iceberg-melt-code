function berg_densities = estimate_iceberg_density(berg_flipped,berg_elev,berg_elevMAD,density_z,density,density_profile,wetdensity_profile)
if berg_flipped == 0 %berg is a fragment or overturned so assume the full firn is saturated
    j=1;
    Hberg(j) = (rho_sw/(rho_sw-rho_i))*berg_elev;
    Hberg_err(j) = abs(Hberg(j)).*sqrt((rho_sw_err/rho_sw)^2 + ((1.4826*berg_elevMAD)/berg_elev)^2 + ((rho_sw_err^2 + rho_i_err^2)/(rho_sw-rho_i)^2));
    wet_ref = find(density_z<=density.eightthir,1,'last');
    if Hberg(j) > density_z(wet_ref)
        rho_prof = [wetdensity_profile(1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j))+1)];
        rho_profrange(1,:) = [wetdensity_profile(2,1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
        rho_profrange(2,:) = [wetdensity_profile(3,1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
    else
        rho_prof = [wetdensity_profile(1:ceil(Hberg(j))+1)];
        rho_profrange(1,:) = [wetdensity_profile(2,1:ceil(Hberg(j))+1)];
        rho_profrange(2,:) = [wetdensity_profile(3,1:ceil(Hberg(j))+1)];
    end
    rho_f(j) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
    rho_f_err(j) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]),[],'omitnan');
    clear rho_prof*;
    %iterate
    while j
        Hberg(j+1) = (rho_sw/(rho_sw-rho_f(j)))*berg_elev;
        Hberg_err(j+1) = abs(Hberg(j+1)).*sqrt(((abs(rho_sw/(rho_sw-rho_f(j)))*sqrt((rho_sw_err/rho_sw)^2 + (sqrt(rho_sw_err^2+rho_f_err(j)^2)/(rho_sw-rho_f(j)))^2))/(rho_sw/(rho_sw-rho_f(j))))^2 + ((1.4826*berg_elevMAD)/berg_elev)^2);
        if Hberg(j+1) > density_z(wet_ref)
            rho_prof = [wetdensity_profile(1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
            rho_profrange(1,:) = [wetdensity_profile(2,1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
            rho_profrange(2,:) = [wetdensity_profile(3,1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
        else
            rho_prof = [wetdensity_profile(1:ceil(Hberg(j+1))+1)];
            rho_profrange(1,:) = [wetdensity_profile(2,1:ceil(Hberg(j+1))+1)];
            rho_profrange(2,:) = [wetdensity_profile(3,1:ceil(Hberg(j+1))+1)];
        end
        rho_f(j+1) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
        rho_f_err(j+1) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]),[],'omitnan');
        %test for convergence
        if abs(rho_f(j+1)-rho_f(j)) < 0.25*rho_f_err(j+1)
            berg_densities(1) = rho_f(j+1);
            berg_densities(2) = nanmean(rho_profrange(1,:));
            berg_densities(3) = nanmean(rho_profrange(2,:));
            clear Hberg rho_f* dry_ref wet_ref rho_prof*;
            break
        else
            j = j+1; clear rho_prof*;
        end
    end
else %berg is upright so only wet the firn below the waterline
    j=1;
    Hberg(j) = (rho_sw/(rho_sw-rho_i))*berg_elev;
    Hberg_err(j) = abs(Hberg(j)).*sqrt((rho_sw_err/rho_sw)^2 + ((1.4826*berg_elevMAD)/berg_elev)^2 + ((rho_sw_err^2 + rho_i_err^2)/(rho_sw-rho_i)^2));
    dry_ref = find(density_z<=berg_elev,1,'last'); wet_ref = find(density_z<=density.eightthir,1,'last');
    if Hberg(j) > density_z(wet_ref)
        rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j))+1)];
        rho_profrange(1,:) = [mindensity_profile(1:dry_ref) wetdensity_profile(2,dry_ref+1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
        rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) wetdensity_profile(3,dry_ref+1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j))+1)];
    else
        rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:ceil(Hberg(j))+1)];
        rho_profrange(1,:) = [mindensity_profile(1:dry_ref) wetdensity_profile(2,dry_ref+1:ceil(Hberg(j))+1)];
        rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) wetdensity_profile(3,dry_ref+1:ceil(Hberg(j))+1)];
    end
    rho_f(j) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
    rho_f_err(j) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]),[],'omitnan');
    clear rho_prof*;
    %iterate
    while j
        Hberg(j+1) = (rho_sw/(rho_sw-rho_f(j)))*berg_elev;
        Hberg_err(j+1) = abs(Hberg(j+1)).*sqrt(((abs(rho_sw/(rho_sw-rho_f(j)))*sqrt((rho_sw_err/rho_sw)^2 + (sqrt(rho_sw_err^2+rho_f_err(j)^2)/(rho_sw-rho_f(j)))^2))/(rho_sw/(rho_sw-rho_f(j))))^2 + ((1.4826*berg_elevMAD)/berg_elev)^2);
        if Hberg(j+1) > density_z(wet_ref)
            rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:wet_ref) density_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
            rho_profrange(1,:) = [mindensity_profile(1:dry_ref) wetdensity_profile(2,dry_ref+1:wet_ref) mindensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
            rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) wetdensity_profile(3,dry_ref+1:wet_ref) maxdensity_profile(wet_ref+1:ceil(Hberg(j+1))+1)];
        else
            rho_prof = [density_profile(1:dry_ref) wetdensity_profile(dry_ref+1:ceil(Hberg(j+1))+1)];
            rho_profrange(1,:) = [mindensity_profile(1:dry_ref) wetdensity_profile(2,dry_ref+1:ceil(Hberg(j+1))+1)];
            rho_profrange(2,:) = [maxdensity_profile(1:dry_ref) wetdensity_profile(3,dry_ref+1:ceil(Hberg(j+1))+1)];
        end
        rho_f(j+1) =  nanmean(rho_prof); %rho_f(j) = rho_i+(f.b*(rho_i-f.a)*(exp(-Hberg(j)/f.b)-1))/Hberg(j); %commented equation is average of the exponential equation for the dry density profile
        rho_f_err(j+1) = max(abs([nanmean(rho_profrange(1,:))-rho_f(j) nanmean(rho_profrange(2,:))-rho_f(j)]),[],'omitnan');
        %test for convergence
        if abs(rho_f(j+1)-rho_f(j)) < 0.25*rho_f_err(j+1)
            berg_densities(1) = rho_f(j+1);
            berg_densities(2) = nanmean(rho_profrange(1,:));
            berg_densities(3) = nanmean(rho_profrange(2,:));
            clear Hberg rho_f* dry_ref wet_ref rho_prof*;
            break
        else
            j = j+1; clear rho_prof*;
        end
    end
end


end