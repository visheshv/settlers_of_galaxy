function idx = find_closest_momentum_star(star_positions_target,star_velocities_target,x0,n,star_ID)

star_id_settled = star_ID(star_ID~=0); % IDs for settled stars

if length(star_id_settled)~=0
star_positions_target(star_id_settled,1:3)= repmat([inf inf inf],length(star_id_settled),1);
end

% Identify 10 closest stars
r0 =  x0(1:3);
v0= x0(4:6);
idx = knnsearch(star_positions_target,r0','K',10);          

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

% stm0 = zeros(1,36); stm0(1,1)=1; stm0(1,8)=1; stm0(1,15)=1; stm0(1,22)=1; stm0(1,29)=1; stm0(1,36)=1;
% ic = horzcat(x0',stm0)';
% xdot = dynamics(0,ic);
% acc_r0 = xdot(4:6);
% velocity_query=v0+acc_r0 * tof + r0 / tof;
% velocity_database= star_velocities_target + (1/tof) * star_positions_target;
% 
% idx = knnsearch(velocity_database,velocity_query');          % Closest star

