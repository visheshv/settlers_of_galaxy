function [idx, t_departure,t_arrival,is_bad_solution,x0_departure] = find_n_best_stars_strategy2(x0,n,star_ID, t_min_departure,query_type,current_star_ID,J_merit,delV_max,delV_used,x,y,z,vx,vy,vz,star_data)

% Input
% starID: current settled Star IDs,
% t: minimum departure
% n: number of target stars
% query_type: 1 for mothership, 2 for settler ship

% current_star_ID: find the ID for the current star to be departed ( !only for settler ships! )

% Output
% t_min_departure
% tof
% Target IDs for the 'n' best stars
% is_bad_solution: equals 0 if legitimately good solution found
% x0 state at updated departure conditions


kpc2km = 30856775814671900;
myr = 1e6*31557600;
kms2kpcpmyr = myr/kpc2km;
kpcpmyr2kms = kpc2km / myr;

counter=0;
angle_thresh_plane=5; % deg
angle_thresh_anomaly=10; % deg
r_thresh=5; % kpc

poly_n=[1.64965819444119e-08,-1.70149835795047e-06,6.54110993481119e-05,-0.00110272848583409,0.00589360129400973,0.0467622383963570,-0.668482140993974,3.15728858911639,-2.33996016100402];

max_check=200;
J_query_store=zeros(max_check,10); % t_dep
stm0 = zeros(1,36); stm0(1,1)=1; stm0(1,8)=1; stm0(1,15)=1; stm0(1,22)=1; stm0(1,29)=1; stm0(1,36)=1; % STM0 for all propagation requirements

for i_counter=1:max_check
    
    % required
    J_merit_temp =0;
    delV_max_temp =0;
    delV_used_temp =0;
    star_ID_temp = star_ID;
    
    i_tof_margin = 20/0.5; % myr
    i_departure = randi([t_min_departure/0.5+1 179],1,1);
    i_arrival = randi([i_departure+1 min(i_departure+i_tof_margin,181)],1,1);
    t_departure = (i_departure-1) * 0.5;
    t_arrival = (i_arrival-1) *0.5;
    tof = t_arrival - t_departure;
    
    %% Assign initial positions at t_departure and not t_departure_min
    if query_type == 1
        
        if (t_departure-t_min_departure)~=0
            tspan = 0:0.1:(t_departure-t_min_departure);
            ic = horzcat(x0',stm0);
            opts = odeset('AbsTol',1e-5);
            [~,states]=ode45(@dynamics,tspan,ic,opts);
            r0 = states(end,1:3)';                     %state
            v0=  states(end,4:6)';
            x0_departure = [r0;v0];
        else
            r0=x0(1:3);
            v0=x0(4:6);
            x0_departure = [r0;v0];
        end
        
    elseif query_type ==2
        
        r0= [x(current_star_ID+1,i_departure),y(current_star_ID+1,i_departure),z(current_star_ID+1,i_departure)]';
        v0= [vx(current_star_ID+1,i_departure),vy(current_star_ID+1,i_departure),vz(current_star_ID+1,i_departure)]';
        x0_departure = [r0;v0];
%         tof_fit=polyval(poly_n,norm(r0));
    end
    
    %% Find stars to target
        
    star_positions_target=[x(2:end,i_arrival),y(2:end,i_arrival),z(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof
    star_velocities_target=[vx(2:end,i_arrival),vy(2:end,i_arrival),vz(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof
    
    star_id_settled = star_ID(star_ID~=0); % IDs for settled stars
    
    % Plane check based on orbit normal direction
    star_momenta=cross(star_positions_target,star_velocities_target);
    star_momenta_dir=star_momenta./vecnorm(star_momenta')';
    star_positions_target_dir=star_positions_target./vecnorm(star_positions_target')';
    star_positions_target_radius=vecnorm(star_positions_target')';
    sc_momentum_vec=cross(r0(1:3),v0(1:3));
    sc_momentum_dir=sc_momentum_vec./vecnorm(sc_momentum_vec')';   % direction of rotation
    sc_momentum_dir=repmat(sc_momentum_dir',1e5,1);
    sc_position_dir=repmat(r0(1:3)'/norm(r0(1:3)),1e5,1);
    sc_position_norm=repmat(norm(r0(1:3)),1e5,1);
    
       
    % Check angle from the initial position
    sc_angle =( 180/pi)* (t_arrival-t_departure) * norm(v0(1:3))/norm(r0(1:3));  % angle (deg)
    
    % Check the radial position 
    idx_plane = find ((abs(star_positions_target_radius-sc_position_norm)<r_thresh) & (acosd(dot(star_momenta_dir,sc_momentum_dir,2))<angle_thresh_plane) & (acosd(dot(star_positions_target_dir,sc_position_dir,2))< (sc_angle+angle_thresh_anomaly)) & (acosd(dot(star_positions_target_dir,sc_position_dir,2))> (sc_angle-angle_thresh_anomaly))); 
    
    idx_plane = setdiff(idx_plane,star_id_settled); % Remove the departing star from search
    
    num_temp=randi([1 length(idx_plane)],1,1);
    i_temp=idx_plane(num_temp);
    
    x_t = [x(i_temp+1,i_arrival) y(i_temp+1,i_arrival) z(i_temp+1,i_arrival) vx(i_temp+1,i_arrival) vy(i_temp+1,i_arrival) vz(i_temp+1,i_arrival)]'; % Target states
    vt = x_t(4:6);rt = x_t(1:3);
    
    
    %% Solve for the dVs
    v0_guess=v0;
    [~,v0,vf]= shooting_solver(r0,rt,tof,v0_guess); % Run the algorithm
    
    delv_transfer= norm(v0-v0_guess');
    delv1= v0-v0_guess';
    delv_rendezvous= norm(-vf + vt');
    delv2= (-vf + vt'); % rendezvous delv
    violation_dv=0;
    
    if query_type ==1
        dV_max_ms= 200; % kmps
        delV_max_temp = delV_max + dV_max_ms; % Add the mothership dv_limits (km/s)
        delV_used_temp = delV_used + delv_transfer * kpcpmyr2kms; % Add the mothership dv_limits (km/s)
        if delv_transfer * kpcpmyr2kms > dV_max_ms
            violation_dv =1;
        end
        
    elseif query_type ==2
        dV_max_ss=175; %kmps
        delV_max_temp = min(delV_max + 2 * dV_max_ss, delV_max +400)  ; % Add the mothership dv_limits (km/s)
        delV_used_temp = delV_used + delv_transfer * kpcpmyr2kms; % Add the mothership dv_limits (km/s)
        if delv_transfer * kpcpmyr2kms > dV_max_ss || delv_rendezvous * kpcpmyr2kms > dV_max_ss
            violation_dv =1;
        end
    end
    
    %% Check  whether it is worth selecting the star with respect to the
    % merit function increase
    
    star_ID_temp = [star_ID(star_ID~=0) ; i_temp];
    error_J_term_temp = J_N_r_theta(star_ID_temp,star_data);
    J_merit_temp = error_J_term_temp * delV_max_temp / delV_used_temp;
    err_term = J_merit * (delV_used /delV_max);
    
    if delV_used_temp < delV_max_temp && violation_dv==0 && error_J_term_temp >= err_term
        J_query_store(i_counter,:) = [i_temp, t_departure,t_arrival, J_merit_temp, x0_departure'];
    else
        J_query_store(i_counter,:) = [i_temp, t_departure,t_arrival, -inf, x0_departure'];
    end
    
end

J_query_store=unique(J_query_store,'rows');
J_query_store=sortrows(J_query_store,4,'descend');

idx= J_query_store(1:n,1);
t_departure= J_query_store(1:n,2);
t_arrival= J_query_store(1:n,3);
x0_departure= J_query_store(1:n,5:10)';

is_bad_solution=sum(J_query_store(1:n,4)<0); % equals 0 if good solution found


end