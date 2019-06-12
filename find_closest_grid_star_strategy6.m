function idx = find_closest_grid_star_strategy6(x0,n,star_ID,x,y,z,t_arrival_temp,star_grid_id,R_vec)

%Inputs
% t_departure_temp: departure time from Sol (myr)
% t_arrival_temp: arrival time from Sol (myr)
% star_grid_id: list of id of grid center stars
% x,y,z,vx,vy,vz: Position and velocity database of all stars
% R_vec: radius database

%Outputs
% idx: star to reach with the fastest tof

% Define initial states
r0 =  x0(1:3);
v0 =  x0(4:6);

star_id_settled = [star_ID(star_ID~=0); 0]; % IDs for settled stars, include sol
star_target_ids=setdiff(star_grid_id,star_id_settled);  % remove settled stars from the target list

% Identify the current position grid center radius
r_upper=min(norm(r0)-mod(norm(r0),1)+3,32);

% select below 5 only if you are close to five
% if norm(r0) >5 && norm(r0) <5.5
%     r_lower=max(norm(r0)-mod(norm(r0),1)-2,2);  % Permission to look for targets below granted
% else
    r_lower=max(norm(r0)-mod(norm(r0),1)-2,5);
% end

i_arrival_temp=t_arrival_temp/0.5 +1;
star_positions_target=[x(1:end,i_arrival_temp),y(1:end,i_arrival_temp),z(1:end,i_arrival_temp)]; % Except sun, all position values for stars at t=tof

star_positions_target(R_vec>r_upper | R_vec<r_lower,:) =  repmat([inf inf inf],sum(R_vec>r_upper | R_vec<r_lower),1);

% Eliminate grid centers which lie beyond the radius band
rel_pos= star_positions_target(star_target_ids+1,:)-repmat(r0',length(star_target_ids),1);
normvec= vecnorm((rel_pos)')';
rel_pos= rel_pos./normvec;

v_n= repmat(v0'/norm(v0),length(star_target_ids),1);

angles = acosd(dot(rel_pos,v_n,2));
angles(isnan(angles))=180;
sign_data=angle_sign(rel_pos,v_n,repmat(r0',length(star_target_ids),1));

angles= angles .* sign_data;

% Identify the star center closest to the current star and having a
% separation angle between (+-90deg)
i_angles=find((angles>-90)&(angles<90));

if isempty(i_angles)
    
    idx=[];
    
else
    
    [sorted_relpos,i_position]=sort(normvec(i_angles));
    
    if isnan(sorted_relpos)
        idx=[];
    else
        i_star_candidate= star_target_ids(i_angles(i_position(1)));
        idx=i_star_candidate;
        
    end
    
end
end