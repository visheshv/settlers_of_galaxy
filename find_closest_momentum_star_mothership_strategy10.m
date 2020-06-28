function N_pool = find_closest_momentum_star_mothership_strategy10(kms2kpcpmyr,R_vec,phi_vec,omega_vec,i_vec,n_vec,v_vec)

% Inputs
% star_positions_target: star positions at arrival (kpc)for star ID 1 to 1e5
% star_velocities_target: star velocities at arrival (kpc)for star ID 1 to 1e5
% x0: initial state 
% n: number of star targets queried
% star_ID: matrix of occupied star IDs
% x,y,z: position database for all stars
% inc: inclination database for star ID 1 to 1e5
% mothership_id: <-
% mothership_controls: matrix containing guidance for mothership {id, t1(1st impulse),t2,t3 (3rd impulse), t4 (max time), max dv1, max dv2, max dv3, max angle1, max angle2, max angle3}
% num_impulses: number of impulses already exhausted by mothership before selection of the current star
% kms2kpcpmyr: conversion to kpc/myr

% Output
% N_pool: candidate starID list


N_pool=[];
iter_count=0;

while(1)
    
    iter_count=iter_count+1;
    % Choose random dt1, dv1, dt2
    dt1= 0.1+rand*10;
    dt2= 2+rand*10;
    dv_mag= 190*rand;
    angle_dv= -40+80*rand;
    
    % Define mothership state at t=dt1
    x0= R_vec(1)*(cosd(n_vec(1)*dt1+phi_vec(1))*cosd(omega_vec(1)) - sind(n_vec(1)*dt1+phi_vec(1))*cosd(i_vec(1))*sind(omega_vec(1)));
    y0= R_vec(1)*(cosd(n_vec(1)*dt1+phi_vec(1))*sind(omega_vec(1)) + sind(n_vec(1)*dt1+phi_vec(1))*cosd(i_vec(1))*cosd(omega_vec(1)));
    z0= R_vec(1)*(sind(n_vec(1)*dt1+phi_vec(1))*sind(i_vec(1)));
    vx0= v_vec(1)*(-sind(n_vec(1)*dt1+phi_vec(1))*cosd(omega_vec(1)) - cosd(n_vec(1)*dt1+phi_vec(1))*cosd(i_vec(1))*sind(omega_vec(1)));
    vy0= v_vec(1)*(-sind(n_vec(1)*dt1+phi_vec(1))*sind(omega_vec(1)) + cosd(n_vec(1)*dt1+phi_vec(1))*cosd(i_vec(1))*cosd(omega_vec(1)));
    vz0= v_vec(1)*(cosd(n_vec(1)*dt1+phi_vec(1))*sind(i_vec(1)));
    
    r0=[x0 y0 z0]';
    v0=[vx0 vy0 vz0]';
    
    dv_vec= dv_mag * [cosd(angle_dv) sind(angle_dv) 0; -sind(angle_dv) cosd(angle_dv) 0; 0 0 1]*v0/norm(v0);  
    x0=[r0;dv_vec+v0];
    
    tspan=dt1:0.1:(dt1+dt2);
    stm0 = zeros(1,36); stm0(1,1)=1; stm0(1,8)=1; stm0(1,15)=1; stm0(1,22)=1; stm0(1,29)=1; stm0(1,36)=1;
    ic = horzcat(x0',stm0); opts = odeset('AbsTol',1e-5);
    [t_hist,states]=ode45(@dynamics,tspan,ic,opts); r_hist=states(:,1:3);v_hist=states(:,4:6);
    
    min_sep=0.1;% minimum separation allowed for the target stars from any of the settled stars
    
    [store_close_stars]=multi_intercept_check2(t_hist,r_hist,v_hist,R_vec,phi_vec,omega_vec,i_vec,v_vec,n_vec,min_sep);

    idx=unique(store_close_stars(:,1));
    
    % inclination of targets
    if idx~=0 & ~isempty(idx)
        inc=180-i_vec(idx+1);
        idx=idx(idx~=0 & inc <3);
        N_pool=[N_pool; idx];
    end
    
    if length(N_pool)>10 | iter_count>10000
        break
    end
 
end

% Check:
% plot3(settled_star_pos_temp(:,1),settled_star_pos_temp(:,2),settled_star_pos_temp(:,3),'o');
% hold on; plot3(x(idx+1,181),y(idx+1,181),z(idx+1,181),'.'); hold on; plot3(x(idx+1,181),y(idx+1,181),z(idx+1,181),'*')
% plot3(star_positions_target(idx,1),star_positions_target(idx,2),star_positions_target(idx,3),'.')
% hold on; plot3(star_positions_target(idx,1),star_positions_target(idx,2),star_positions_target(idx,3),'*')
