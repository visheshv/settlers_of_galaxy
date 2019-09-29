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
global J_merit delV_max delV_used
global R_min R_max
global settled_star_pos

%Constants
kpc2km = 30856775814671900;
myr = 1e6*31557600;
kms2kpcpmyr = myr/kpc2km;
kpcpmyr2kms = kpc2km / myr;
dtr = pi/180;
R_min = 2;
R_max = 32;

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
i_arrival=t_arrival/0.5+1;                  %column index corresponding to t = 6myr
tof= t_arrival-t_departure;                 %myr
x0 = [x(1,i_departure) y(1,i_departure) z(1,i_departure) vx(1,i_departure) vy(1,i_departure) vz(1,i_departure)]';  % Sol position at t0
r0 = x0(1:3);
v0_guess = x0(4:6); % Guess initial condition for solver

num_impulses_ms=0;
sum_impulses_ms=0;
store_results =zeros(1,16);
delv_store=[];
state_store=[];
settlement_tree_sp=[];
settled_star_pos = [];

J_merit =0;
delV_max = 0;
delV_used = 0;
J_merit_temp =0;
delV_max_temp =0;
delV_used_temp =0;

tic

constraints_met=0; % check whether delta V constraints are met after finding a solution

while constraints_met==0
    
    % Sol-less data
    star_positions_target=[x(2:end,i_arrival),y(2:end,i_arrival),z(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof
    star_velocities_target=[vx(2:end,i_arrival),vy(2:end,i_arrival),vz(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof
    
    idx = find_closest_momentum_star_strategy1(star_positions_target,star_velocities_target,x0,1,star_ID,x,y,z);
    x_t = [x(idx+1,i_arrival) y(idx+1,i_arrival) z(idx+1,i_arrival) vx(idx+1,i_arrival) vy(idx+1,i_arrival) vz(idx+1,i_arrival)]'; % Target states
    
    vt = x_t(4:6);
    rt = x_t(1:3);
    
    [states,v0,vf]= shooting_solver(r0,rt,tof,v0_guess);
    
    store_results = [x_t(1:6)' states(1,4:6) states(end, 4:6) x_t(4:6)'  tof];
    
    delv_transfer= norm(v0-v0_guess');
    delv1= v0-v0_guess';
    delv_rendezvous= vecnorm((store_results(:,10:12)-store_results(:,13:15))')';
    delv2= (-store_results(:,10:12)+store_results(:,13:15)); % rendezvous delv
    
    % Check  whether it is worth selecting the star with respect to the
    % merit function increase
    delV_max_temp = delV_max + 300; % Add the mothership dv_limits (km/s)
    delV_used_temp = delV_used + delv_transfer * kpcpmyr2kms + delv_rendezvous * kpcpmyr2kms; % Add the mothership dv_limits (km/s)
    star_ID_temp = [star_ID(star_ID~=0) ; idx];
    error_J_term_temp = J_N_r_theta(star_ID_temp,star_data);
    J_merit_temp = error_J_term_temp * delV_max_temp / delV_used_temp;
    
    if delv_transfer * kpcpmyr2kms > 200 || delv_rendezvous * kpcpmyr2kms > 300  % If Mothership initial impulse exceeds  200 km/s or setlling pod impulse exceeds 300 km/s
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
        
        % Update mothership impulses
        num_impulses_ms=num_impulses_ms+1;
        sum_impulses_ms=sum_impulses_ms+delv_transfer;
        
        if num_impulses_ms>3 || sum_impulses_ms * kpcpmyr2kms > 500 % Dont use more than three impulses or more than 500 kmps dV
            % Generation 0 mothership star captures and update the star
            % settled data base
            %             star_database_ms = find_solution_intercepts_ms(t_departure,r0,v0_guess,idx_prev);
            break
        elseif num_impulses_ms==1
            % Update settled star ids
            star_ID(1,num_impulses_ms) = idx; % 1st settlement from Mothership
            delV_max = delV_max + 500; % Add the mothership dv_limits (km/s)
            error_J_term = J_N_r_theta(star_ID(star_ID~=0),star_data);
            J_merit = error_J_term * delV_max / delV_used;
        else
            if norm(rt)-norm(r0)<0
               break
            end
            % Update settled star ids
            star_ID(1,num_impulses_ms) = idx; % 1st settlement from Mothership
            error_J_term = J_N_r_theta(star_ID(star_ID~=0),star_data);
            J_merit = error_J_term * delV_max / delV_used;
        end
        
        delv_store = [delv_store; delv1 delv_transfer delv2 delv_rendezvous t_departure t_arrival];
        state_store = [state_store; states];
        settlement_tree_sp =[settlement_tree_sp; -1 num_impulses_ms idx t_arrival delv2* kpcpmyr2kms];
        
        % Update mothership states
        tof = 2.5;                                %myr
        t_departure=t_arrival+1.5;                %myr since maneuver has to happen >1 yr after SP maneuver
        i_departure=t_departure/0.5+1;            %column number query
        
        % Propagate to t_departure
        stm0 = zeros(1,36); stm0(1,1)=1; stm0(1,8)=1; stm0(1,15)=1; stm0(1,22)=1; stm0(1,29)=1; stm0(1,36)=1;
        tspan = 0:0.1:1.5;
        ic = horzcat(states(end,1:6),stm0);
        opts = odeset('AbsTol',1e-5);
        [~,states]=ode45(@dynamics,tspan,ic,opts);
        
        x0 = states(end,1:6)';                     %state
        
        t_arrival= t_departure+tof;               %myr
        i_arrival=t_arrival/0.5+1;                        % column index corresponding to t = 5myr
        tof= t_arrival-t_departure;               %myr
        
        r0 =  x0(1:3);
        v0_guess = x0(4:6); % Guess initial condition for solver
        idx_prev = idx;
        
    end

end

%% add lines for extra motherships1,2 , NOTE: not adding the mothership matrix now
settlement_tree_sp =[settlement_tree_sp; -2 1 53794 14.5 4.12468957441322,166.542907556157,102.042781906309];
star_ID(1,4)=53794;

toc;
toc-tic
% Mothership line
settlement_tree_ms(1,:)=[-1, 0, length(star_ID(star_ID~=0)), length(star_ID(star_ID~=0)), delv_store(:,end-1)', reshape(delv_store(:,1:3)',1,[]) * kpcpmyr2kms];
fprintf(fileID,[repmat('%0.0f,',1,4) repmat('%0.12f,',1,12)],settlement_tree_ms);
% Settlement Pod line
fprintf(fileID,['\n' repmat('%0.0f,',1, 3) repmat('%0.12f,',1,4)],settlement_tree_sp');

%% fast ship line 
star_ID(1,10:11) = [11916 19131]; % Solutions available from J_store

%!! RUN the fast_ship_transfer file for FS -11 
settlement_tree_fs(1:2,:)=[-11,11916,0,24.5,45.9127551644999,803.175220136587,-14.3951226842251,319.705806225252,-614.358394676057,31.3572788494068;
                           -12,19131,0,6.5,-1316.70466616146,362.474264571612,9.85548837500831,1382.67201066288,-333.627000067822,-24.4211899177846];  % row 1 and row3 of initial solution
% FS1 to theta =-90 deg : -11,696,0,25,kpcpmyr2kms*[-0.607263816523397,-0.365907897190692,-0.0446800262498599,0.408446251323942,0.635341048639990,-0.00469600342441436
fprintf(fileID,['\n' repmat('%0.0f,',1,2) repmat('%0.12f,',1,8)],settlement_tree_fs');


%% Greedy strategy for settlement
% For departure time of 8 (6+2) myr, identify the three closest stars
% positionally and make a Gen 1 set of settler ships

tic
gen = 1; % 1st generation of settler ships to be seen
settlement_tree_ss=[];

while(gen<22)

% remove the stars that are already occupied for generating the search
% space
if gen==7
    disp('Check for gen7 issues')    
end
disp(['Current gen:' num2str(gen+1)])

idx_search = star_data(2:end,1);       % Complete search space for stars except sol
star_id_settled = star_ID(star_ID~=0); % IDs for settled stars

% idx_search = setdiff(idx_search,star_id_settled); % Remove settled star IDs from search space

% Set of upto 3 ^ k stars at Gen k, find 3 closest star for each settled
% star of Gen k such that upto 3 ^ (k+1) stars selected. 
stars_gen_k = star_ID(gen,:);
stars_gen_k = stars_gen_k(stars_gen_k~=0); % find the latest generation settled stars
%r_temp = star_positions_target;

% Make the positions of the settled stars high enough to be rejected by
% the distance query
% for k = 1 : length(star_id_settled)
%         id_temp = find(star_id_settled(k)==star_data(2:end,1));
%         r_temp(id_temp+1,1:3) = [+inf,+inf,+inf];
% end

% Search star pairs for all stars in the gen k
for  j =1 : length(stars_gen_k)
        
    tof = 2.5; % starting with small guess
    ID_jk = stars_gen_k(j); % star ID of the jth element of the kth gen
    
    if gen==1
        
        if j<=(length(stars_gen_k)-2)
            t_departure=settlement_tree_sp(j,4)+2;    % 2 yr for settling a star system
        else
            t_departure=settlement_tree_fs(j+2-length(stars_gen_k),4)+2;    % 2 yr for settling a star system
        end
        
        t_arrival= t_departure+tof;               % myr
        i_arrival = t_arrival/0.5 +1;             % position history column corresponding to the arrival time
        i_departure = t_departure/0.5+1;          % position history column corresponding to the departure time
        star_positions_target=[x(2:end,i_arrival),y(2:end,i_arrival),z(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof
        
    else
        
        indx_current_star = find(settlement_tree_ss(:,2) == ID_jk);
        t_departure=settlement_tree_ss(indx_current_star(1),5)+2; % 2 yr
        t_arrival= t_departure+tof;                 %myr
        i_arrival = t_arrival/0.5 +1;             % position history column corresponding to the arrival time
        i_departure = t_departure/0.5+1;          % position history column corresponding to the departure time
        
        if t_arrival>89.5
            continue
        end
        
        star_positions_target=[x(2:end,i_arrival),y(2:end,i_arrival),z(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof           
    end
    
    
    idx_jk = find(ID_jk == star_data(:,1)); % find the row number for the ID in the star database
    r0_jk = [x(idx_jk,i_departure) y(idx_jk,i_departure) z(idx_jk,i_departure)];   % Position of the star ID of the jth element of the kth gen at departure time
    v0_jk = [vx(idx_jk,i_departure) vy(idx_jk,i_departure) vz(idx_jk,i_departure)];   % Position of the star ID of the jth element of the kth gen at departure time
    x0 = [r0_jk';v0_jk'];
    idx_vec = find_closest_momentum_star_strategy1(star_positions_target,star_velocities_target,x0,3,star_ID,x,y,z) ; % 3 closest stars settlerships
        
     for l =1:length(idx_vec)  % Find all the closest 1/2/3 stars with respect to jth settled star of kth generation
        
        t_arrival_temp=t_arrival;
        t_departure_temp=t_departure;
        i_arrival_temp = t_arrival_temp/0.5 +1;             % position history column corresponding to the arrival time
        i_departure_temp = t_departure_temp/0.5+1;          % position history column corresponding to the departure time
                
        constraints_met =0;
        num_constraint_violation =0;
        
        J_merit_temp =0;
        delV_max_temp =0;
        delV_used_temp =0;
        
        while constraints_met == 0
            idx = idx_vec(l);
            x_t = [x(idx+1,i_arrival_temp) y(idx+1,i_arrival_temp) z(idx+1,i_arrival_temp) vx(idx+1,i_arrival_temp) vy(idx+1,i_arrival_temp) vz(idx+1,i_arrival_temp)]'; % Target states
            
            % Solve the shooting problem for all star pairs and remove the pairs which
            % are problematic in terms of transfer deltaV, rendezvous deltaV,  tof =inf
            tof = t_arrival_temp-t_departure_temp;
            r0 = r0_jk';
            rt = x_t(1:3);
            vt = x_t(4:6);
            v0_guess = [vx(idx_jk,i_departure_temp) vy(idx_jk,i_departure_temp) vz(idx_jk,i_departure_temp)]'; % guess value provided as the star velocity itself
            
            [states,v0,vf]= shooting_solver(r0,rt,tof,v0_guess);
            
            dv1= v0-v0_guess'; % transfer impulse
            delv_transfer = norm(dv1);
            dv2= -vf+vt';      % rendezvous impulse
            delv_rendezvous = norm(dv2);
            
            % Check  whether it is worth selecting the star with respect to the
            % merit function increase
            delV_max_temp = delV_max + 400; % Add the mothership dv_limits (km/s)
            delV_used_temp = delV_used + delv_transfer * kpcpmyr2kms + delv_rendezvous * kpcpmyr2kms; % Add the mothership dv_limits (km/s)
            star_ID_temp = [star_ID(star_ID~=0) ; idx];
            error_J_term_temp = J_N_r_theta(star_ID_temp,star_data);
            J_merit_temp = error_J_term_temp * delV_max_temp / delV_used_temp;
            
            if norm(rt)-norm(r0)<0
               break
            end
            
            if norm(dv1) * kpcpmyr2kms > 175 || norm(dv2) * kpcpmyr2kms > 175 % || J_merit_temp -J_merit <0 % no solution
%                 disp(['no solution at tof (myr):' num2str(tof)])
                tof =tof+0.5;
                t_arrival_temp= t_arrival_temp+0.5;                 %myr
                i_arrival_temp=t_arrival_temp/0.5+1;                          %column index corresponding to t = 6myr
                num_constraint_violation = num_constraint_violation+1;
                
                if num_constraint_violation> 40 || t_arrival_temp>89.5
                    break
                end
                
            else 
                
                % Successful transfer
                delV_used = delV_used + delv_transfer * kpcpmyr2kms + delv_rendezvous * kpcpmyr2kms; % km/s
                delV_max = delV_max + 400; % km/s
                error_J_term = J_N_r_theta(star_ID(star_ID~=0),star_data);
                J_merit = error_J_term * delV_max / delV_used;
                
                star_ID(gen+1,(j-1)*3+l)=star_data(idx+1,1);
                % Settlement Pod line
                % Write the actions in the txt solution file. 
                settlement_tree_ss=[settlement_tree_ss; ID_jk idx 2 t_departure_temp t_arrival_temp dv1* kpcpmyr2kms dv2* kpcpmyr2kms];
                fprintf(fileID,['\n' repmat('%0.0f,', 1, 3) repmat('%0.12f,', 1, 8) ],settlement_tree_ss(end,:));
                break
                
            end
            
        end
    end
end

toc
% Update the departure and arrival times, generation gen
gen = gen+1;     % next generation of settler ships to be seen

% if gen>7
% save 4jundata
% end

toc-tic
end


close all
% plot star positions versus generation
for i = 1:9
    for j=1:length(star_ID(star_ID~=0))
        i_f=90/0.5+1;
        if star_ID(i,j)==0
        else
%         temp = find(star_ID(i,j) == star_data(:,1)); % find the row number for the ID in the star database
%         idx_ij =  star_data(temp,1);
        idx_ij = star_ID(i,j)+1;
        x_t = [x(idx_ij,i_f) y(idx_ij,i_f) z(idx_ij,i_f) vx(idx_ij,i_f) vy(idx_ij,i_f) vz(idx_ij,i_f)]'; % Target states
        figure(1)
        hold on
        plot3(x_t(1),x_t(2),x_t(3),'*')
        axis equal
        xlim([-32 32])
        ylim([-32 32])
        end
    end
%     pause(5)    
end
 hold on; plot3(0,0,0,'+')
