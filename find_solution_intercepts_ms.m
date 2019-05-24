function star_database_ms = find_solution_intercepts_ms(t_departure,r0,v0_guess,idx)

% Find out possible solutions for star interceptions after the mothership
% has exhausted its last impulse
global R_vec phi_vec omega_vec i_vec n_vec


star_database_ms=[];
x0 = [r0;v0_guess];  % Sol position at t0
stm0 = zeros(1,36); stm0(1,1)=1; stm0(1,8)=1; stm0(1,15)=1; stm0(1,22)=1; stm0(1,29)=1; stm0(1,36)=1;
t_remaining= 90 - t_departure;
tspan = 0:0.05:t_remaining;
ic = horzcat(x0',stm0);
[t,states]=ode45(@dynamics,tspan,ic);

parfor i = 1:length(t)
    
    ms_pos=states(i,1:3);
    time= t(i) + t_departure;
    pos_stars=zeros(length(R_vec),3);
    
    for j=1:length(R_vec)
        
    x= R_vec(j)*(cosd(n_vec(j)*time+phi_vec(j))*cosd(omega_vec(j)) - sind(n_vec(j)*time+phi_vec(j))*cosd(i_vec(j))*sind(omega_vec(j)));
    y= R_vec(j)*(cosd(n_vec(j)*time+phi_vec(j))*sind(omega_vec(j)) + sind(n_vec(j)*time+phi_vec(j))*cosd(i_vec(j))*cosd(omega_vec(j)));
    z= R_vec(j)*(sind(n_vec(j)*time+phi_vec(j))*sind(i_vec(j)));
    pos_stars(j,:)= [x y z];
    
    end
    
    pos_stars(idx+1,:)=[inf,inf, inf];  % Make sure that initial star is not taken as an interception
    
    idx_min = knnsearch(pos_stars,ms_pos);
    min_dist= norm(pos_stars(idx_min,:)-ms_pos);
    
    if min_dist<1e-4 && norm(ms_pos)>2-0.01 && norm(ms_pos)<32-0.01
       star_database_ms=[star_database_ms; time idx_min+1];
    end
    
end