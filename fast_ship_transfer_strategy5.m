function [id,t_arrival,delv1,delv2]=fast_ship_transfer_strategy5(t_departure,r_search_min,r_search_max,theta_min,theta_max,x,y,z,vx,vy,vz,i_vec)

%Inputs
% t_departure: departure time from Sol (myr)
% r_search_min: min radius of search(kpc)
% r_search_max: max radius of search(kpc)
% theta_search: min theta of search(deg) at the time of arrival
% theta_search: max radius of search(deg) at the time of arrival

%Outputs
% id: star to reach with the fastest tof
% t_arrival: arrival time (myr)
% delv1: Transfer dV (kmps)
% delv2: Rendezvous dV (kmps)

%Constants
kpc2km = 30856775814671900; myr = 1e6*31557600;  kpcpmyr2kms = kpc2km / myr;

% load all the star positions greater than 15 kpc with tof =2.5 and
% solve for the transfers with tof++ until fast ship transfers are
% acceptable.
i_departure=t_departure/0.5+1;      %column number query
t_arrival= t_departure+2.5;         %myr
i_arrival=t_arrival/0.5+1;          %column index corresponding to t_arrival

star_positions_target=[x(2:end,i_arrival),y(2:end,i_arrival),z(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof

inc=180-i_vec;

% Filter stars based on position
r_stars=(vecnorm(star_positions_target'))';
idx_fs_targets= find((r_stars>r_search_min) & (r_stars<r_search_max) & inc<5);


tic
store=[];

% Iterate for all stars under idx_fs_targets

for j=1:length(idx_fs_targets)
    
    i_arrival=t_arrival/0.5+1;          %column index corresponding to t_arrival
    tof= t_arrival-t_departure;
    x0 = [x(1,i_departure) y(1,i_departure) z(1,i_departure) vx(1,i_departure) vy(1,i_departure) vz(1,i_departure)]';  % Sol position at t0
    r0 = x0(1:3); v0_guess = x0(4:6); % Guess initial condition for solver
    
    while(1)
        
        idx = idx_fs_targets(j);
        x_t = [x(idx+1,i_arrival) y(idx+1,i_arrival) z(idx+1,i_arrival) vx(idx+1,i_arrival) vy(idx+1,i_arrival) vz(idx+1,i_arrival)]'; % Target states
        
        theta_star= atan2d(y(idx+1,i_arrival),x(idx+1,i_arrival));
        
        if theta_star>theta_max || theta_star<theta_min
            break
        end
                   
        vt = x_t(4:6);
        rt = x_t(1:3);
        
        [states,v0,vf]= shooting_solver(r0,rt,tof,v0_guess);
        
        delv_transfer= norm(v0-v0_guess');
        delv1= v0-v0_guess';
        delv_rendezvous= norm(-vf+vt');
        delv2= (-vf+vt'); % rendezvous delv
        
        r_min=min(vecnorm((states(:,1:3))')');
        
        if (delv_transfer * kpcpmyr2kms + delv_rendezvous * kpcpmyr2kms) > 1500 || r_min < 2 % If Fast ship total impulse exceeds  1500 km/s
            tof =tof+0.5;
            t_arrival= t_departure+tof;                         %myr
            i_arrival=t_arrival/0.5+1;                          %column index corresponding to t = 6myr            
            if t_arrival>89.5
                break
            end            
        else
            store=[store; idx tof (delv_transfer * kpcpmyr2kms + delv_rendezvous * kpcpmyr2kms) r_stars(idx)  theta_star delv1*kpcpmyr2kms delv2*kpcpmyr2kms t_arrival];
            break
        end
        
    end
end

if ~isempty(store)
    store=sortrows(store,2);    
    id=store(1,1);              % starID
    t_arrival=store(1,end);     % myr
    delv1 = store(1,5:7);       % kpc/myr
    delv2 = store(1,8:10);      % kpc/myr
    
else
    id =0; % incase there are no solutions
    disp('Could not find fast ship solutions')    
end






