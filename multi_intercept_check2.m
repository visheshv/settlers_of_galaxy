function [store_close_stars]=multi_intercept_check2(t_hist,r_hist,v_hist,R_vec,phi_vec,omega_vec,i_vec,v_vec,n_vec,min_sep)

% Inputs 
% t_hist: time history (n x 1) [sec]
% r_hist: position history (n x 3) [kpc]
% v_hist: vel history (n x 3) [kpc/myr]
% R_vec,phi_vec,omega_vec,i_vec,v_vec,n_vec: 100k star database radius (kpc), anomaly (deg), rate(deg/s), inclination(deg), speed(kpc/myr), ang frequency

% Step 1 Filter stars that have crossed the mothership trajectory

store_close_stars=zeros(length(t_hist),2);

r0_ms=r_hist(1,:)';
rf_ms=r_hist(end,:)';
zcap=[0 0 1]';
n_cap=cross(zcap,(rf_ms-r0_ms)/norm(rf_ms-r0_ms));
n_cap=n_cap/norm(n_cap);

% figure(1)
% plot3(r_hist(:,1),r_hist(:,2),r_hist(:,3),'.')
% hold on

store_crossing_stars=zeros(length(R_vec),1);

for i=1:length(R_vec)
        
        t0= t_hist(1);tf=t_hist(end);
        x0= R_vec(i)*(cosd(n_vec(i)*t0+phi_vec(i))*cosd(omega_vec(i)) - sind(n_vec(i)*t0+phi_vec(i))*cosd(i_vec(i))*sind(omega_vec(i)));
        y0= R_vec(i)*(cosd(n_vec(i)*t0+phi_vec(i))*sind(omega_vec(i)) + sind(n_vec(i)*t0+phi_vec(i))*cosd(i_vec(i))*cosd(omega_vec(i)));
        z0= R_vec(i)*(sind(n_vec(i)*t0+phi_vec(i))*sind(i_vec(i)));
        xf= R_vec(i)*(cosd(n_vec(i)*tf+phi_vec(i))*cosd(omega_vec(i)) - sind(n_vec(i)*tf+phi_vec(i))*cosd(i_vec(i))*sind(omega_vec(i)));
        yf= R_vec(i)*(cosd(n_vec(i)*tf+phi_vec(i))*sind(omega_vec(i)) + sind(n_vec(i)*tf+phi_vec(i))*cosd(i_vec(i))*cosd(omega_vec(i)));
        zf= R_vec(i)*(sind(n_vec(i)*tf+phi_vec(i))*sind(i_vec(i)));    
        
        r0=[x0 y0 z0]';
        rf=[xf yf zf]';
        
        ang1= acosd(dot((r0-r0_ms)/norm(r0-r0_ms),n_cap));
        ang2= acosd(dot((rf-r0_ms)/norm(rf-r0_ms),n_cap));
        
        if norm(r0)>norm(r0_ms)*0.98 && norm(r0)<norm(rf_ms)*1.02 && ((ang1 <90 && ang2 >90) ||  (ang1 >90 && ang2 <90))
            store_crossing_stars(i)=1;
        end
        
       
end

idx_cross_stars=find(store_crossing_stars==1);

for j=1:length(t_hist)

    r=r_hist(j,:);
    t=t_hist(j);
    v=v_hist(j,:);
    
    star_positions_target=zeros(1e5+1,3);
    star_vel_target=zeros(1e5+1,3);
    
    for i=1:length(R_vec)
    
        if i==1
            star_positions_target(i,:)=[inf inf inf];
            star_vel_target(i,:)=[inf inf inf];
        else
            
            x= R_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*cosd(omega_vec(i)) - sind(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*sind(omega_vec(i)));
            y= R_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*sind(omega_vec(i)) + sind(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*cosd(omega_vec(i)));
            z= R_vec(i)*(sind(n_vec(i)*t+phi_vec(i))*sind(i_vec(i)));
            
            star_positions_target(i,:)=[x y z];

        end
        
    end
    
    idx = knnsearch(star_positions_target,r,'K',1);     
    dist_min=norm(star_positions_target(idx,:)- r);
    
    if dist_min<min_sep
        store_close_stars(j,:)= idx-1;  % -1 to identify the star ID
    end
    
end

    

