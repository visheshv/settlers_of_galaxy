function idx = find_n_best_stars(star_positions_target,star_velocities_target,x0,n,star_ID, t_min_departure)

% Input
% starID: current settled Star IDs, 
% t: minimum departure
% n:n target stars

% Output
% t_min_departure
% tof
% Target IDs for the 'n' best stars 

global J_merit

star_id_settled = star_ID(star_ID~=0); % IDs for settled stars

if length(star_id_settled)~=0
star_positions_target(star_id_settled,1:3)= repmat([inf inf inf],length(star_id_settled),1);
end

% Identify 10 closest stars
r0 =  x0(1:3);
v0= x0(4:6);
idx = knnsearch(star_positions_target,r0','K',100);          

% Find the closest of them in the direction of the current velocity
rel_pos= star_positions_target(idx,:)-repmat(r0',length(idx),1);
normvec= vecnorm((rel_pos)')';
rel_pos= rel_pos./normvec;
v_n= repmat(v0'/norm(v0),length(idx),1);

angles = acosd(dot(rel_pos,v_n,2));
% [~,i_angle]=min(angles);
[~,i_angles]=sort(angles);

idx=idx(i_angles);
idx=idx(1:n);