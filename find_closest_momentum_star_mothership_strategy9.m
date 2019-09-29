function idx_selection = find_closest_momentum_star_mothership_strategy9(x0,n,star_ID,x,y,z, mothership_id, mothership_controls, num_impulses,kms2kpcpmyr,R_vec,phi_vec,omega_vec,i_vec,n_vec,v_vec,t_departure)

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
% idx: target starID


% Define initial states
r0 =  x0(1:3);
v0 =  x0(4:6);

min_sep=1;% minimum separation allowed for the target stars from any of the settled stars

star_id_settled = star_ID(star_ID~=0); % IDs for settled stars

if ~isempty(star_id_settled)
    settled_star_pos_temp=[x(star_id_settled+1,181) y(star_id_settled+1,181) z(star_id_settled+1,181)];       % temporary matrix to store settled star positions at tf =90 myr
end

% Add selection constraints based on mothership controls
dv_allowed = mothership_controls(-mothership_id,6+num_impulses)*kms2kpcpmyr * 0.9; % Reduce by a factor of 0.9, the max dV to ensure selection
max_angle_dv = mothership_controls(-mothership_id,9+num_impulses);
time_allowed = mothership_controls(-mothership_id,3+num_impulses);
dv_vec= dv_allowed * [cosd(max_angle_dv) sind(max_angle_dv) 0; -sind(max_angle_dv) cosd(max_angle_dv) 0; 0 0 1]*v0/norm(v0);  

% theta_max= sign(max_angle_dv) * acosd(dot((dv_vec+v0)/norm(dv_vec+v0),v0/norm(v0)));% allowed maximum angle to search targets with respect to velocity
 
% check for closest stars
x0=[r0;dv_vec+v0];
tspan=0:0.5:time_allowed;
stm0 = zeros(1,36); stm0(1,1)=1; stm0(1,8)=1; stm0(1,15)=1; stm0(1,22)=1; stm0(1,29)=1; stm0(1,36)=1;
ic = horzcat(x0',stm0); opts = odeset('AbsTol',1e-5);
[t_hist,states]=ode45(@dynamics,tspan,ic,opts); r_hist=states(:,1:3);v_hist=states(:,4:6);

t_hist=t_hist+t_departure; % Add initial time 

[store_close_stars]=multi_intercept_check1(t_hist,r_hist,v_hist,R_vec,phi_vec,omega_vec,i_vec,v_vec,n_vec);
idx=sort(unique(store_close_stars(:,1)));     


% Eliminate angles corresponding to targets which are not at a min separation from settled stars at tf
if ~isempty(star_id_settled)                                        
    for j=1:length(idx)
        r_j = [x(idx(j)+1,181),y(idx(j)+1,181),z(idx(j)+1,181)];
        idx_temp= knnsearch(settled_star_pos_temp,r_j,'K',1);
        dist_min=norm(settled_star_pos_temp(idx_temp,:)-r_j);
        if dist_min < min_sep
            idx(j)=0;
        end
    end
end

% inclination of targets
if idx~=0 | ~isempty(idx)
    inc=180-i_vec(idx+1);
    idx=idx(idx~=0 & inc <3);
        if idx~=0 | ~isempty(idx)
            idx_selection = idx(1);   % Choose any index
        else
            idx_selection=0;
        end
else
    idx_selection=0;
end
 
end

% Check:
% plot3(settled_star_pos_temp(:,1),settled_star_pos_temp(:,2),settled_star_pos_temp(:,3),'o');
% hold on; plot3(x(idx+1,181),y(idx+1,181),z(idx+1,181),'.'); hold on; plot3(x(idx+1,181),y(idx+1,181),z(idx+1,181),'*')
% plot3(star_positions_target(idx,1),star_positions_target(idx,2),star_positions_target(idx,3),'.')
% hold on; plot3(star_positions_target(idx,1),star_positions_target(idx,2),star_positions_target(idx,3),'*')
