function idx = find_closest_momentum_star_strategy1(star_positions_target,star_velocities_target,x0,n,star_ID)

global J_merit

% Define initial states
r0 =  x0(1:3);
v0= x0(4:6);

gen_jump=0;% kpc
% remove stars under or outside the sphere of radius norm(r0)
% condition=vecnorm((star_positions_target)')' >  norm(r0)+gen_jump;
% star_positions_target(condition,:)= repmat([inf inf inf], sum(condition),1 );

star_id_settled = star_ID(star_ID~=0); % IDs for settled stars

if length(star_id_settled)~=0
star_positions_target(star_id_settled,1:3)= repmat([inf inf inf],length(star_id_settled),1);
end

% Identify 'n' closest stars
idx = knnsearch(star_positions_target,r0','K',1000);          
% idx_low = knnsearch(star_positions_target,r0','K',300);
% idx = setdiff(idx_up,idx_low);

rel_pos= star_positions_target(idx,:)-repmat(r0',length(idx),1);
normvec= vecnorm((rel_pos)')';
rel_pos= rel_pos./normvec;
v_n= repmat(v0'/norm(v0),length(idx),1);

angles = acosd(dot(rel_pos,v_n,2));

% Find the highest frequency of them in the direction of the current velocity
% hist=histogram(normvec);[~,i_r_min]=max(hist.Values);
% r_min=hist.BinEdges(i_r_min);
r_min=2;

inc_thresh=2; % inclination range

% plane normal
h=cross(star_positions_target(idx,:),star_velocities_target(idx,:));
h=h./vecnorm((h)')';
inc=acosd(dot(h,repmat([0 0 -1],length(idx),1),2));

angles(normvec<r_min)=repmat(179, sum(normvec<r_min),1 ); % remove the very close stars from contention by making the angle values to be 179
angles(inc>inc_thresh)=repmat(179, sum(inc>inc_thresh),1); % remove the very close stars from contention by making the angle values to be 179

[temp,i_angles]=sort(angles);

if temp(1)==179
   disp('break') 
end

idx=idx(i_angles);
idx=idx(1:n);
