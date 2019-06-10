function idx = find_closest_momentum_star_ss_strategy5(star_positions_target,star_velocities_target,x0,n,star_ID,x,y,z,i_vec,min_sep,r_max,r_min,min_search_angle)

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

if ~isempty(star_id_settled)
    settled_star_pos_temp=[x(star_id_settled+1,181) y(star_id_settled+1,181) z(star_id_settled+1,181)];       % temporary matrix to store settled star positions at tf =90 myr
    star_positions_target(star_id_settled,1:3)= repmat([inf inf inf],length(star_id_settled),1);
end

% Identify 10,000 closest stars
idx = knnsearch(star_positions_target,r0','K',10000);

rel_pos= star_positions_target(idx,:)-repmat(r0',length(idx),1);
normvec= vecnorm((rel_pos)')';
rel_pos= rel_pos./normvec;
rel_vel= star_velocities_target(idx,:)-repmat(v0',length(idx),1);
normvel= vecnorm((rel_vel)')';

v_n= repmat(v0'/norm(v0),length(idx),1);

angles = acosd(dot(rel_pos,v_n,2));
sign_data=angle_sign(rel_pos,v_n,repmat(r0',length(idx),1));

% CHECK!
if norm(r0) >=16  % i.e. r is between 15 and 20 kpc
    
    r_min = 1;
    r_max = r_min + 3;
    min_sep=3;
    
end

if ~isempty(star_id_settled)
    for j=1:length(idx)
        r_j = [x(idx(j)+1,181),y(idx(j)+1,181),z(idx(j)+1,181)];
        idx_temp= knnsearch(settled_star_pos_temp,r_j,'K',1);
        dist_min=norm(settled_star_pos_temp(idx_temp,:)-r_j);
        if dist_min < min_sep
            angles(j) = 179;
        end
    end
end

v_max=0.5;    % velocity threshold

inc_thresh=3; % inclination range acceptable (deg)

% plane normal
inc=180-i_vec(idx);
    
angles(normvec<r_min)=179;
angles(normvec>r_max)=179;
angles(normvel>v_max)=179;
angles(inc>inc_thresh)=179;

angles_temp =  angles .* sign_data;

% find the angle such that min separation between target stars is as
% per the grid size at the radius band

% separation_angle= atand((norm(r0) * pi/16)/ r_min) ; % grid size %% could be also written as  (norm(r0) * pi/16)/ r_min for better calculation
% separation_angle= (180/pi)*((norm(r0) * pi/100)/ r_min) ; % grid size %% could be also written as  (norm(r0) * pi/16)/ r_min for better calculation

% theta_dvdr=atand(dvdr);
% theta_star=max([theta_dvdr,separation_angle]);

% Choose first star to be close to -5 deg from the velocity vector
[~,ind_min]=min(abs(angles_temp-min_search_angle));
idx1=idx(ind_min);
idx=idx1;

% [tempA,i_angles]=sort(angles_temp);
% 
% % Selection strategy: Find one angle as close to radial outwards direction
% % and find the other two angles such that they maximize the separation
% % between the three angles and ensure they are atleast separated by a
% % threshold angle. If not, skip 3 star selection and go for either 1/2
%     
% if separation_angle <= 45
%     
%     i_1= find((tempA> -45) & (tempA< 0));  
%     
%     if length(i_1)>=2
%         if abs(tempA(i_1(1))-tempA(i_1(end)))> 0.8 * separation_angle
%             idx2=idx(i_angles(i_1(1)));
%             idx3=idx(i_angles(i_1(end)));
%             idx=[idx1;idx2;idx3];
%         else
%             idx3=idx(i_angles(i_1(end)));
%             idx=[idx1;idx3];
%         end
%         
%     elseif length(i_1)==1
%         
%         idx2=idx(i_angles(i_1(end)));
%         idx=[idx1;idx2];
%         
%     else
%         idx=idx1;
%     end
%     
% elseif separation_angle > 45 && separation_angle <=90
%     
%     i_1= find((tempA> -90-separation_angle) & (tempA< -90+separation_angle));  
%     
%     if length(i_1)>=2
%         if abs(tempA(i_1(1))-tempA(i_1(end)))> 0.8 * separation_angle  && idx(i_angles(i_1(1)))~=idx1
%             idx2=idx(i_angles(i_1(1)));
%             idx3=idx(i_angles(i_1(end)));
%             idx=[idx1;idx2;idx3];
%         else
%             idx3=idx(i_angles(i_1(end)));
%             idx=[idx1;idx3];
%         end        
%     elseif length(i_1)==1 && idx(i_angles(i_1(end)))~=idx1        
%         idx2=idx(i_angles(i_1(end)));
%         idx=[idx1;idx2];        
%     else
%         idx=idx1;
%     end        
% else
%     
% disp(['Impossible separation angle for SS(deg):' num2str(separation_angle)])
%     
% end

end



% Check:
% plot3(settled_star_pos_temp(:,1),settled_star_pos_temp(:,2),settled_star_pos_temp(:,3),'o');
% hold on; plot3(x(idx+1,181),y(idx+1,181),z(idx+1,181),'.'); hold on; plot3(x(idx+1,181),y(idx+1,181),z(idx+1,181),'*')
% plot3(star_positions_target(idx,1),star_positions_target(idx,2),star_positions_target(idx,3),'.')
% hold on; plot3(star_positions_target(idx,1),star_positions_target(idx,2),star_positions_target(idx,3),'*')
