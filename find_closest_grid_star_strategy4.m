function idx = find_closest_grid_star_strategy4(star_positions_target,star_velocities_target,x0,n,star_ID,star_grid_id,x,y,z,vx,vy,vz,i_arrival,i_departure)

% Define initial states
r0 =  x0(1:3);
v0 = x0(4:6);
v0_guess = v0;

list_all_star_ids= 1:1e5;
star_id_na=setdiff(list_all_star_ids,star_grid_id);
star_id_settled = star_ID(star_ID~=0); % IDs for settled stars

star_id_na=union(star_id_na,star_id_settled);

star_target_ids=setdiff(star_grid_id,star_id_settled);

% Identify 'n' closest stars
dv_store=zeros(length(star_target_ids),2);

parfor i = 1:length(star_target_ids) 
    
    rt=[x(star_target_ids(i)+1,i_arrival),y(star_target_ids(i)+1,i_arrival),z(star_target_ids(i)+1,i_arrival)];
    vt=[vx(star_target_ids(i)+1,i_arrival),vy(star_target_ids(i)+1,i_arrival),vz(star_target_ids(i)+1,i_arrival)];
    
    tof=(i_arrival-i_departure)*0.5;
    
    [~,v0,vf]= shooting_solver(r0,rt,tof,v0_guess);
    
    delv_transfer= norm(v0-v0_guess');
    delv_rendezvous= norm(vf-vt);
    
    if delv_transfer > 175 
    dv_store(i,:)=[star_target_ids(i) delv_transfer+delv_rendezvous];    
    
end

dv_store=unique(dv_store,'rows');
dv_store=sortrows(dv_store,2,'ascend');

if length(dv_store)>=3
    idx = dv_store(1:n,1);
else
    disp('No solutions found')
end


