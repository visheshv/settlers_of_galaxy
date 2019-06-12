function idx = find_closest_momentum_star_ss_strategy6(x0,n,star_ID,x,y,z,i_vec,min_sep,r_max,r_min,min_search_angle,star_grid_id,t_arrival_temp)

% Inputs
% star_positions_target: star positions at arrival (kpc)for star ID 1 to 1e5
% star_velocities_target: star velocities at arrival (kpc)for star ID 1 to 1e5
% x0: initial position
% n: number of star targets queried
% star_ID: matrix of occupied star IDs
% x,y,z: position database for all stars
% i_vec: inclination database for star ID 1 to 1e5
% kms2kpcpmyr: conversion to kpc/myr
% min_sep: minimum separation distance allowed for the target stars from any of the settled stars (kpc)
% min_search_angle: angle away from the velocity vector to search nearby stars (deg)

% Output
% idx: target starID

% Define initial states
r0 =  x0(1:3);
v0= x0(4:6);

% constants, use them to calculate dv/dr
kpc2km = 30856775814671900; myr = 1e6*31557600;
k = [0.00287729 0.0023821 -0.0010625 0.000198502 -1.88428e-05 9.70521e-07 -2.70559e-08 3.7516e-10 -1.94316e-12];
k0 = k(1);k1=k(2);k2=k(3);k3=k(4);k4=k(5);k5=k(6);k6=k(7);k7=k(8);k8=k(9);
rnorm=norm(r0);dvdr=-(myr*(8*k8*rnorm^7 + 7*k7*rnorm^6 + 6*k6*rnorm^5 + 5*k5*rnorm^4 + 4*k4*rnorm^3 + 3*k3*rnorm^2 + 2*k2*rnorm + k1))/(kpc2km*(k8*rnorm^8 + k7*rnorm^7 + k6*rnorm^6 + k5*rnorm^5 + k4*rnorm^4 + k3*rnorm^3 + k2*rnorm^2 + k1*rnorm + k0)^2);

star_id_settled = star_ID(star_ID~=0); % IDs for settled stars

i_arrival_temp=t_arrival_temp/0.5 +1;
star_positions_target=[x(1:end,i_arrival_temp),y(1:end,i_arrival_temp),z(1:end,i_arrival_temp)]; % Except sun, all position values for stars at t=tof

if ~isempty(star_id_settled)
        star_positions_target(star_id_settled+1,1:3)= repmat([inf inf inf],length(star_id_settled),1);
end

% Remove the stars of the non-star database ids 
all_id=1:1e5; 
rem_id=setdiff(all_id,star_grid_id);
star_positions_target(rem_id+1,1:3)= repmat([inf inf inf],length(rem_id),1);
star_positions_target(1,1:3)= repmat([inf inf inf],1,1);

init_num=500;
while_count=0;
min_search_angle_temp= min_search_angle;

while_count_max=25;

while(while_count<while_count_max)

while_count=while_count+1;
    
% Identify 10,000 closest stars
idx = knnsearch(star_positions_target,r0','K',init_num);
rel_pos= star_positions_target(idx,:)-repmat(r0',length(idx),1);
normvec= vecnorm((rel_pos)')';
rel_pos= rel_pos./normvec;

v_n= repmat(v0'/norm(v0),length(idx),1);

angles = acosd(dot(rel_pos,v_n,2));
sign_data=angle_sign(rel_pos,v_n,repmat(r0',length(idx),1));

inc_thresh=10; % inclination range acceptable (deg)

% plane normal
inc=180-i_vec(idx);

if norm(r0)>17
    angles(vecnorm(star_positions_target(idx,:)')'<17)=179;    
end

if norm(r0)>5.5
    angles(vecnorm(star_positions_target(idx,:)')'<5)=179;
end

%%


angles(normvec<r_min)=179;
angles(normvec>r_max)=179;
angles(inc>inc_thresh)=179;

angles_temp =  angles .* sign_data;

% Choose first star to be close to -5 deg from the velocity vector
[~,ind_min]=min(abs(angles_temp-min_search_angle_temp));

if isempty(ind_min)
    
    if mod(while_count,12)==1
        
    r_min=r_min+1;
    r_max = r_max+3;
    init_num= init_num+1000;
    idx=[];
    min_search_angle_temp=min_search_angle;
    
    else
       
        min_search_angle_temp=min_search_angle - 15 * floor(while_count*0.5) * (-1)^while_count;
        if (while_count==while_count_max)
           idx=[];
        end
        
    end
    
    
else
    
    idx1=idx(ind_min);
    idx=idx1-1;
    break
end

end



% Check:
% plot3(settled_star_pos_temp(:,1),settled_star_pos_temp(:,2),settled_star_pos_temp(:,3),'o');
% hold on; plot3(x(idx+1,181),y(idx+1,181),z(idx+1,181),'.'); hold on; plot3(x(idx+1,181),y(idx+1,181),z(idx+1,181),'*')
% plot3(star_positions_target(idx,1),star_positions_target(idx,2),star_positions_target(idx,3),'.')
% hold on; plot3(star_positions_target(idx,1),star_positions_target(idx,2),star_positions_target(idx,3),'*')
