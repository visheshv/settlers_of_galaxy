% Analyse fast ship transfer to various higher R for tof, and polar angle spread
clear all 
close all

load star_snapshots
load star_data

%Constants
kpc2km = 30856775814671900; myr = 1e6*31557600; kms2kpcpmyr = myr/kpc2km; kpcpmyr2kms = kpc2km / myr; dtr = pi/180; R_min = 2; R_max = 32;

% load the sol position at earliest departure time

num_impulses_ms=0; sum_impulses_ms=0; store_results = zeros(1,16);
settlement_tree_sp=[]; settled_star_pos = [];

% load all the star positions greater than 15 kpc with tof =2.5 and
% solve for the transfers with tof++ until fast ship transfers are
% acceptable.
t_departure=0;                      %myr according to CB's analysis for min deltaV
i_departure=t_departure/0.5+1;      %column number query
t_arrival= 2.5;                     %myr
i_arrival=t_arrival/0.5+1;          %column index corresponding to t = 6myr
star_positions_target=[x(2:end,i_arrival),y(2:end,i_arrival),z(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof
star_velocities_target=[vx(2:end,i_arrival),vy(2:end,i_arrival),vz(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof

h=cross(star_positions_target,star_velocities_target);
h=h./vecnorm((h)')';
inc=acosd(dot(h,repmat([0 0 -1],length(h),1),2));

r_stars=(vecnorm(star_positions_target'))';
theta_stars= atan2d(star_positions_target(:,2),star_positions_target(:,1));
idx_fs_targets= find((r_stars>10) & (r_stars<10.3) & (theta_stars<5) & (theta_stars >3) & inc<5);


J_merit =0; delV_max = 0; delV_used = 0; J_merit_temp =0; delV_max_temp =0; delV_used_temp =0;

tic

store=[];

for j=1:length(idx_fs_targets)

% Sol-less data
t_departure=0;                      %myr according to CB's analysis for min deltaV
i_departure=t_departure/0.5+1;      %column number query
t_arrival= 2.5;                     %myr
i_arrival=t_arrival/0.5+1;          %column index corresponding to t = 6myr
tof= t_arrival-t_departure;
x0 = [x(1,i_departure) y(1,i_departure) z(1,i_departure) vx(1,i_departure) vy(1,i_departure) vz(1,i_departure)]';  % Sol position at t0
r0 = x0(1:3); v0_guess = x0(4:6); % Guess initial condition for solver


while(1)
idx = idx_fs_targets(j);
x_t = [x(idx+1,i_arrival) y(idx+1,i_arrival) z(idx+1,i_arrival) vx(idx+1,i_arrival) vy(idx+1,i_arrival) vz(idx+1,i_arrival)]'; % Target states

vt = x_t(4:6);
rt = x_t(1:3);

[states,v0,vf]= shooting_solver(r0,rt,tof,v0_guess);

delv_transfer= norm(v0-v0_guess');
delv1= v0-v0_guess';
delv_rendezvous= norm(-vf+vt');
delv2= (-vf+vt'); % rendezvous delv

r_min=min(vecnorm((states(:,1:3))')');

if (delv_transfer * kpcpmyr2kms) >200 || (delv_rendezvous * kpcpmyr2kms) > 300 || r_min < 2 % If Mothership initial impulse exceeds  200 km/s or setlling pod impulse exceeds 300 km/s
    disp(['no solution at tof (myr):' num2str(tof)])
    tof =tof+0.5;
    t_arrival= t_departure+tof;                 %myr
    i_arrival=t_arrival/0.5+1;                          %column index corresponding to t = 6myr
    
    if t_arrival>89.5
        break
    end
else
    
    % Successful transfer
    delV_used = delV_used + delv_transfer * kpcpmyr2kms + delv_rendezvous * kpcpmyr2kms; % km/s
    delV_max = delV_max + 300; % km/s
    
    if norm(rt)-norm(r0)<0
        break
    end
    
    % Update settled star ids
    break
end

end

store=[store; idx tof (delv_transfer * kpcpmyr2kms + delv_rendezvous * kpcpmyr2kms) r_stars(idx)  theta_stars(idx) delv1*kpcpmyr2kms delv2*kpcpmyr2kms];

end

% store=sortrows(store,2);






