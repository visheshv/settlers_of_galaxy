% Implement a simple strategy to populate stars only with one mothership

% 1. Find closest star of MS for a 0.5 mn tof for which the rendezvous deltaV < 300 km/s
% 2. do loop
%     For each settled star, go forward 2 mn yrs and find three closest stars for each Settler ship, such that initial deltaV transfer is <175 kmps and rendezvous deltaV is <175 kmps
%     while t<90 myr
% 3. Find the settlement tree

clc
clear all
close all

%% Load data
load star_snapshots 
load star_data
 
global star_data  R_vec phi_vec omega_vec i_vec theta_f 
global v_vec n_vec

%Constants
kpc2km = 30856775814671900;
myr = 1e6*31557600;
kms2kpcpmyr = myr/kpc2km;
kpcpmyr2kms = kpc2km / myr;
dtr = pi/180;

% Update star database variables
R_vec = data(:,2); % kpc (kiloparsecs) 
phi_vec = data(:,5); % phi (deg)
omega_vec = data(:,4); % omega (deg)
i_vec = data(:,3); % i (deg)
theta_f= data(:,6); % final polar angle (deg)
k = [0.00287729 0.0023821 -0.0010625 0.000198502 -1.88428e-05 9.70521e-07 -2.70559e-08 3.7516e-10 -1.94316e-12];
v_vec = kms2kpcpmyr*([R_vec.^0 R_vec.^1 R_vec.^2 R_vec.^3 R_vec.^4 R_vec.^5 R_vec.^6 R_vec.^7 R_vec.^8]*k').^-1; % Star speeds kpc/Myr
n_vec = (1/dtr)*(v_vec ./ R_vec); % (deg/mYR)

star_data = data; % Load the star IDs into this variable
star_ID   = zeros(100,1e5); % Assuming 100 generations of settlements amd 100000 stars to be potentially covered

fileID = fopen('strategy1.txt','w');
fprintf(fileID,'%s\n','strategy1');

% Outputs
% analysis of dV vs tof. starting from Sol, find and transfer to the closest
% star at each snap

%% Mother ship capture of near star with minimum transfer deltaV
% Initialize sol state

t_departure=2.5;                    %myr according to CB's analysis for min deltaV
i_departure=t_departure/0.5+1;      %column number query
t_arrival= 5;                       %myr
i=t_arrival/0.5+1;                  %column index corresponding to t = 6myr
tof= t_arrival-t_departure;                 %myr
x0 = [x(1,i_departure) y(1,i_departure) z(1,i_departure) vx(1,i_departure) vy(1,i_departure) vz(1,i_departure)]';  % Sol position at t0
r0 = x0(1:3);
v0_guess = x0(4:6); % Guess initial condition for solver

num_impulses_ms=0;
sum_impulses_ms=0;
store_results =zeros(1,16);
delv_store=[];
state_store=[];
tic

constraints_met=0; % check whether delta V constraints are met after finding a solution

while constraints_met==0
    
    % Sol-less data
    star_positions_target=[x(2:end,i),y(2:end,i),z(2:end,i)]; % Except sun, all position values for stars at t=tof
    star_velocities_target=[vx(2:end,i),vy(2:end,i),vz(2:end,i)]; % Except sun, all position values for stars at t=tof
    
    idx = find_closest_momentum_star(star_positions_target,star_velocities_target,x0);
    x_t = [x(idx+1,i) y(idx+1,i) z(idx+1,i) vx(idx+1,i) vy(idx+1,i) vz(idx+1,i)]'; % Target states
    
    vt = x_t(4:6);
    rt = x_t(1:3);
    
    star_ID(1,1) = idx; % 1st settlement from Mothership
    
    [states,v0,vf]= shooting_solver(r0,rt,tof,v0_guess);
    
    store_results = [x_t(1:6)' states(1,4:6) states(end, 4:6) x_t(4:6)'  tof];
    
    delv_transfer= norm(v0-v0_guess');
    delv1= v0-v0_guess';
    delv_rendezvous= vecnorm((store_results(:,10:12)-store_results(:,13:15))')';
    delv2= (-store_results(:,10:12)+store_results(:,13:15)); % rendezvous delv
    
    if delv_transfer * kpcpmyr2kms > 200 || delv_rendezvous * kpcpmyr2kms > 300  % If Mothership initial impulse exceeds  200 km/s or setlling pod impulse exceeds 300 km/s
        disp(['no solution at tof (myr):' num2str(tof)])
        tof =tof+0.5;
        t_arrival= t_departure+tof;                 %myr
        i=t_arrival/0.5+1;                          %column index corresponding to t = 6myr
    else
        
        % Update mothership impulses
        num_impulses_ms=num_impulses_ms+1;
        sum_impulses_ms=sum_impulses_ms+delv_transfer;
        
        if num_impulses_ms>3 || sum_impulses_ms * kpcpmyr2kms > 500 % Dont use more than three impulses or more than 500 kmps dV
            % Generation 0 mothership star captures and update the star
            % settled data base
            star_database_ms = find_solution_intercepts_ms(t_departure,r0,v0_guess,idx_prev);
            break 
        end
        
        delv_store = [delv_store; delv1 delv_transfer delv2 delv_rendezvous];
        state_store = [state_store; states];
        
        % Update mothership states
        t_departure=t_arrival;                    %myr according to CB's analysis for min deltaV
        i_departure=t_departure/0.5+1;            %column number query
        t_arrival= t_departure+tof;                 %myr
        i=t_arrival/0.5+1;                        % column index corresponding to t = 6myr
        tof= t_arrival-t_departure;               %myr
        x0 = x_t;
        r0 =  x_t(1:3);
        v0_guess = x_t(4:6); % Guess initial condition for solver
        idx_prev = idx;
           
    end

end

toc;
toc-tic
% Mothership line
settlement_tree_ms(1,:)=[-1, 0, 1, 1, 2.5, delv_store(end,1:3) * kpcpmyr2kms];
fprintf(fileID,'%d\t%d\t%d\t%d\t%6.5f\t%6.5f\t%6.5f\t%6.5f\t\n',settlement_tree_ms);
% Settlement Pod line
settlement_tree_sp=[-1, 1, star_ID(1,1), 6, delv_store(end,5:7)* kpcpmyr2kms];
fprintf(fileID,'%d\t%d\t%d\t%6.5f\t%6.5f\t%6.5f\t%6.5f\n',settlement_tree_sp);

%% Greedy strategy for settlement
% For departure time of 8 (6+2) myr, identify the three closest stars
% positionally and make a Gen 1 set of settler ships

gen = 1; % 1st generation of settler ships to be seen
t_departure = 8; % myr
t_arrival = 12; % 1st generation tof set as 2 myr;

while(gen<12)

% remove the stars that are already occupied for generating the search
% space

idx_search = star_data(2:end,1);       % Complete search space for stars
star_id_settled = star_ID(star_ID~=0); % IDs for settled stars
idx_search = setdiff(idx_search,star_id_settled); % Remove settled star IDs from search space
i_arrival = t_arrival/0.5 +1;             % position history column corresponding to the arrival time
i_departure = t_departure/0.5+1;          % position history column corresponding to the departure time

star_positions_target=[x(2:end,i_arrival),y(2:end,i_arrival),z(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof

% Set of upto 3 ^ k stars at Gen k, find 3 closest star for each settled
% star of Gen k such that upto 3 ^ (k+1) stars selected. 
stars_gen_k = star_ID(gen,:);
stars_gen_k = stars_gen_k(stars_gen_k~=0); % find the latest generation settled stars
r_temp = star_positions_target;

% Make the positions of the settled stars high enough to be rejected by
% the distance query
for k = 1 : length(star_id_settled)
        id_temp = find(star_id_settled(k)==star_data(2:end,1));
        r_temp(id_temp+1,1:3) = [+inf,+inf,+inf];
end

% Search star pairs for all stars in the gen k
for  j =1 : length(stars_gen_k)
    
    ID_jk = stars_gen_k(j); % star ID of the jth element of the kth gen
    idx_jk = find(ID_jk == star_data(:,1)); % find the row number for the ID in the star database
    r0_jk = [x(idx_jk,i_departure) y(idx_jk,i_departure) z(idx_jk,i_departure)];   % Position of the star ID of the jth element of the kth gen at departure time
    idx_vec = knnsearch(r_temp,r0_jk,'K',3);                                       % 3 Closest stars
    
    for l =1:3
    
        idx = idx_vec(l);
        x_t = [x(idx+1,i_arrival) y(idx+1,i_arrival) z(idx+1,i_arrival) vx(idx+1,i_arrival) vy(idx+1,i_arrival) vz(idx+1,i_arrival)]'; % Target states
        
        % Solve the shooting problem for all star pairs and remove the pairs which
        % are problematic in terms of transfer deltaV, rendezvous deltaV,  tof =inf
        tof = t_arrival-t_departure;
        r0 = r0_jk';
        rt = x_t(1:3);
        vt = x_t(4:6);
        v0_guess = [vx(idx_jk,i_departure) vy(idx_jk,i_departure) vz(idx_jk,i_departure)]'; % guess value provided as the star velocity itself
        
        [states,v0,vf]= shooting_solver(r0,rt,tof,v0_guess);
        
        dv1= v0-v0_guess'; % transfer impulse
        dv2= -vf+vt';      % rendezvous impulse
        
        if norm(dv1) * kpcpmyr2kms > 175 || norm(dv2) * kpcpmyr2kms > 175 % no solution
%             dv1_applied = 175* kms2kpcpmyr *v0/norm(v0);
%             v0_restricted = v0_guess' + dv1_applied;
%             tspan = 0:0.1:1;
%             
%             % make a function whose inputs are (r0,rt,tof,v0) and gives out
%             % the five impulses required by the settlership
%             
%             
%             [~,states]=ode45(@dynamics,tspan,ic);
%             
%             v0_guess = states(end,4:6)';             % Update v0 for the second part solution
%             r0 = states(end,1:3)';                   % 
%             tof = t_arrival-t_departure-1;
            star_ID(gen+1,(j-1)*3+l)=star_data(idx+1);
        else
            star_ID(gen+1,(j-1)*3+l)=star_data(idx+1);
            % Update the settled star database
            r_temp(idx,1:3) = [+inf,+inf,+inf]; % Make the target star out of reach for the j+1th element            
        end
        
    end
    
% Write the actions in the txt solution file. 


end
% Update the departure and arrival times, generation gen
gen = gen+1;     % next generation of settler ships to be seen
t_departure = t_arrival+2; % myr
t_arrival = t_departure+4; % 1st generation tof set as 2 myr;

end

% plot star positions versus generation
for i = 1:9
    for j=1:3^(i-1)
        i_arrival = (12+6 *(i-1))/0.5+1;
        temp = find(star_ID(i,j) == star_data(:,1)); % find the row number for the ID in the star database
        idx_ij =  star_data(temp,1);
        x_t = [x(idx_ij,i_arrival) y(idx_ij,i_arrival) z(idx_ij,i_arrival) vx(idx_ij,i_arrival) vy(idx_ij,i_arrival) vz(idx_ij,i_arrival)]'; % Target states
        figure(1)
        hold on
        plot3(x_t(1),x_t(2),x_t(3),'*')
    end
    
end

