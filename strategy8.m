% Strategy 6: Important points

% Implement a simple strategy to populate three stars with three motherships
% Implement star selection strategy (1,2,3) with 3-impulse shooter for 3.5 kpc stars. Implement the shooter. FS1 used here. Must make sure that radial ships targets are atleast
% Three motherships are made to evolve at various angles to cover max path. One mothership can leave at 10 myr to cover maximum polar angle of around 220 deg.
% Fast ship goes towards the inertial north east  (r~15-20) as fast as possible.

clc
% clearvars -except iter_count
clear all
close all

%% Load data
load star_snapshots
load star_data
load star_target_database_inc_10


%% globals
global star_data  R_vec phi_vec omega_vec i_vec theta_f
global v_vec n_vec
global J_merit delV_max delV_used
global R_min R_max
global x y z vx vy vz

warning('off','all')

%% Constants
kpc2km = 30856775814671900;
myr = 1e6*31557600;
kms2kpcpmyr = myr/kpc2km;
kpcpmyr2kms = kpc2km / myr;
dtr = pi/180;
R_min = 2; R_max = 32;

% For selecting tof, tof polynomial function versus r (kpc)
poly_n=[1.64965819444119e-08,-1.70149835795047e-06,6.54110993481119e-05,-0.00110272848583409,0.00589360129400973,0.0467622383963570,-0.668482140993974,3.15728858911639,-2.33996016100402];

%% Update star database variables
R_vec = data(:,2); phi_vec = data(:,5); omega_vec = data(:,4); i_vec = data(:,3); theta_f= data(:,6); % Units kpc, deg, deg/myr, deg,deg
k = [0.00287729 0.0023821 -0.0010625 0.000198502 -1.88428e-05 9.70521e-07 -2.70559e-08 3.7516e-10 -1.94316e-12];
v_vec = kms2kpcpmyr*([R_vec.^0 R_vec.^1 R_vec.^2 R_vec.^3 R_vec.^4 R_vec.^5 R_vec.^6 R_vec.^7 R_vec.^8]*k').^-1; % Star speeds kpc/myr
n_vec = (1/dtr)*(v_vec ./ R_vec); % (deg/mYR)

star_data = data; % Load the star IDs into this variable
star_ID   = zeros(100,1e5); % Assuming 100 generations of settlements amd 100000 stars to be potentially covered

J_merit =0; delV_max = 0; delV_used = 0;

%% Enter test input vector
vv = [0.00000	10.00000	5.00000	15.00000	100.00000	100.00000	20.00000	20.00000	-30.00000	90.00000	5.00000	10.00000	10.00000	5.00000	100.00000	50.00000	30.00000	0.00000	90.00000	80.00000	5.00000	10.00000	5.00000	5.00000	150.00000	80.00000	50.00000	50.00000	90.00000	120.00000	27  27.1  -90  0    27   27.1  0  90	1.00000	3.00000	1.00000	-5.00000];

itr = '1';


%% Solution format text file
fname = 'strategy8';
fnameext = strcat(fname,'.txt');
figname = strcat(fname, '.fig');
fileID = fopen(fnameext,'w'); fprintf(fileID,'%s','strategy8');

settlement_tree_sp=[];
settlement_tree_ms=[];

%% Control points
% Controls for motherships: id, t_departure (1st impulse),t1 (time period of first segment after departure burn,t2 (myr), t3 (myr), max dv1 (km/s), max dv2 (km/s), max dv3 (km/s), max angle1 (deg), max angle2 (deg), max angle3 (deg)
mothership_controls=[-1, vv(1:10);
    -2, vv(11:20);
    -3, vv(21:30)];

% Control for fast ships
% ArgIn: t_departure: departureSol(myr), r_search_min: min radius of search(kpc), r_search_max: max radius of search(kpc),theta_search: min theta of search(deg) at arrival, theta_search: max radius of search(deg) at the time of arrival
t_departure_fs1=0; r_query_min_fs1=vv(31); r_query_max_fs1=vv(32); theta_query_min_fs1=vv(33); theta_query_max_fs1=vv(34);
t_departure_fs2=0; r_query_min_fs2=vv(35); r_query_max_fs2=vv(36); theta_query_min_fs2=vv(37); theta_query_max_fs2=vv(38);

% Control for separation from existing settled stars(distance) for settler ships
min_sep= vv(39); % kpc
r_max= vv(40);   % kpc
r_min= vv(41);   % kpc
min_search_angle= vv(42); % deg


%% Mother ship capture of near star with minimum transfer deltaV
for mothership_id=-1:-1:-3
    
    t_departure=mothership_controls(-mothership_id,2);     %myr
    
    num_impulses_ms=0; sum_impulses_ms=0;                  % Update the impulse count to zero for each mothership at the initial iteration start
    i_departure=t_departure/0.5+1;                         %column number query
    t_arrival= t_departure+2.5;                            %myr
    i_arrival=t_arrival/0.5+1;                             %column index corresponding to t = 6myr
    tof= t_arrival-t_departure;                            %myr
    
    x0 = [x(1,i_departure) y(1,i_departure) z(1,i_departure) vx(1,i_departure) vy(1,i_departure) vz(1,i_departure)]';  % Sol position at tdeparture
    r0 = x0(1:3); v0_guess = x0(4:6); % Sol velocity at departure, Guess initial condition for solver
    
    
    store_results =zeros(1,16);
    delv_store=[]; state_store=[];
    
    J_merit_temp =0;delV_max_temp =0;delV_used_temp =0;
    
    tic
    constraints_met=0; % check whether delta V constraints are met after finding a solution
    num_violations=0;
    
    while constraints_met==0
        
        star_not_found_tof=0;
        % data of star positions except Sol
        star_positions_target=[x(2:end,i_arrival),y(2:end,i_arrival),z(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof
        star_velocities_target=[vx(2:end,i_arrival),vy(2:end,i_arrival),vz(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof
        
        if num_impulses_ms==3
            break
        end
        
        idx=find_closest_momentum_star_mothership_strategy5(x0,1,star_ID,x,y,z, mothership_id, mothership_controls, num_impulses_ms,kms2kpcpmyr,R_vec,phi_vec,omega_vec,i_vec,n_vec,v_vec,t_departure);
        
        if length(idx)>3
            disp('Problem identifying stars')
            break
        end
        
        if idx~=0
            r0=x0(1:3);v0_guess= x0(4:6);
            
            x_t = [x(idx+1,i_arrival) y(idx+1,i_arrival) z(idx+1,i_arrival) vx(idx+1,i_arrival) vy(idx+1,i_arrival) vz(idx+1,i_arrival)]'; % Target states
            rt = x_t(1:3);vt = x_t(4:6);
            
            tof=t_arrival-t_departure;
            
            [states,v0,vf]= shooting_solver(r0,rt,tof,v0_guess);                       % Run the two point shooting algorithm
            
            delv_transfer= norm(v0-v0_guess');
            delv1= v0-v0_guess';                                                       % transfer dV (kpc/myr)
            delv_rendezvous= norm(-vf+vt');
            delv2= (-vf+vt');                                                          % rendezvous delv (kpc/myr)
        else
            star_not_found_tof=1;
        end
        
        %% Identify successful or failed transfer
        if delv_transfer * kpcpmyr2kms > 200 || delv_rendezvous * kpcpmyr2kms > 300 || star_not_found_tof==1 % If Mothership initial impulse exceeds  200 km/s or setlling pod impulse exceeds 300 km/s
            
            tof =tof+0.5;
            t_arrival= t_departure+tof;                 %myr
            i_arrival=t_arrival/0.5+1;                          %column index corresponding to t = 6myr
            
            if t_arrival>89.5
                break                                           %break the while loop because there is no successful star position
            end
            
            num_violations=num_violations+1;
            
            if num_violations>20
                disp(['No good star found at the desired angle. Mothership,Impulse number,StarID:' num2str(mothership_id) ',' num2str(num_impulses_ms) ',' num2str(idx)])
                break
            end
            
        else
            
            % Successful transfer
            delV_used = delV_used + delv_transfer * kpcpmyr2kms + delv_rendezvous * kpcpmyr2kms; % km/s
            delV_max = delV_max + 300; % km/s % For the settling pods
            
            % Update mothership impulses
            num_impulses_ms=num_impulses_ms+1;
            sum_impulses_ms=sum_impulses_ms+delv_transfer;
            
            if num_impulses_ms>3 || sum_impulses_ms * kpcpmyr2kms > 500 || (t_arrival+1.5)>89.5 % Dont use more than three impulses or more than 500 kmps dV
                disp(['No good star found at the desired angle. Mothership,Impulse number,StarID, totaldV (km/s), t_arrival(myr):' num2str(mothership_id) ',' num2str(num_impulses_ms) ',' num2str(idx) ',' num2str(sum_impulses_ms * kpcpmyr2kms) ',' num2str(t_arrival)])
                break
            elseif num_impulses_ms==1
                star_ID(1,num_impulses_ms+ 3*(-mothership_id-1)) = idx; % 1st settlement from Mothership
                delV_max = delV_max + 500;                              % Add the mothership dv_limits (km/s) only once
            else
                star_ID(1,num_impulses_ms+ 3*(-mothership_id-1)) = idx; % nth settlement from Mothership
            end
            
            % Update merit function terms
            error_J_term = J_N_r_theta(star_ID(star_ID~=0),star_data);
            J_merit = error_J_term * delV_max / delV_used;
            
            delv_store = [delv_store; t_departure delv1* kpcpmyr2kms t_arrival delv2* kpcpmyr2kms mothership_id];
            settlement_tree_sp =[settlement_tree_sp; mothership_id num_impulses_ms idx t_arrival delv2* kpcpmyr2kms];
            
            % Update mothership states
            t_margin=mothership_controls(-mothership_id,2+num_impulses_ms)-t_arrival;     % Identify the time left to do the next burn myr
            t_sp_margin= 1.5;                                                          % Minimum time to next burn
            
            t_departure=t_arrival+max(t_sp_margin,t_margin);                           %myr since maneuver has to happen >1 yr after SP rendezvous maneuver
            
            % Propagate to t_departure which is ahead of the t_arrival by
            % more than 1 myr
            stm0 = zeros(1,36); stm0(1,1)=1; stm0(1,8)=1; stm0(1,15)=1; stm0(1,22)=1; stm0(1,29)=1; stm0(1,36)=1;
            tspan = 0:0.1:max(t_sp_margin,t_margin);
            ic = horzcat(states(end,1:6),stm0);
            opts = odeset('AbsTol',1e-5);
            [~,states]=ode45(@dynamics,tspan,ic,opts);
            
            i_departure=t_departure/0.5+1;                         %column number query
            t_arrival= t_departure+2.5;                            %myr
            i_arrival=t_arrival/0.5+1;                             %column index corresponding to t = 6myr
            
            % Update initial state for the next star capture
            x0 = states(end,1:6)';                     %state
            r0 =  x0(1:3);                             %position
            v0_guess = x0(4:6);                        %Guess initial condition for solver
            idx_prev = idx;
            
        end
    end
    
    %% Mothership line
    length_ms_info=size(delv_store,1);                                                                                  % number of rows of delv_store = number of impulses
    
    if length_ms_info==0
        continue
    else
        dvrow=reshape(delv_store(:,2:4)',1,numel(delv_store(:,2:4)));                                                       % dv's put in a row
        settlement_tree_ms=[mothership_id, 0, length_ms_info, length_ms_info, delv_store(:,1)',dvrow];                      % number of dvs= number of SPs
        fprintf(fileID,['\n' repmat('%d,',1,4) repmat('%0.12f,',1,4*length_ms_info)],settlement_tree_ms(end,:)');
        
        % Settlement Pod line
        fprintf(fileID,['\n' repmat('%0.0f,',1, 3) repmat('%0.12f,',1,4)],settlement_tree_sp(end-length_ms_info+1:end,:)'); % last three rows of the mothership
        
    end
end

toc;
toc-tic;

%% Fastship update solution file
% Fastship line
[id_fs1,t_arrival_fs1,delv1_fs1,delv2_fs1]=fast_ship_transfer_strategy7(t_departure_fs1,r_query_min_fs1,r_query_max_fs1,theta_query_min_fs1,theta_query_max_fs1,x,y,z,vx,vy,vz,i_vec(2:end),star_database);
[id_fs2,t_arrival_fs2,delv1_fs2,delv2_fs2]=fast_ship_transfer_strategy7(t_departure_fs2,r_query_min_fs2,r_query_max_fs2,theta_query_min_fs2,theta_query_max_fs2,x,y,z,vx,vy,vz,i_vec(2:end),star_database);

if id_fs1 ~=0 & id_fs2 ~=0
    settlement_tree_fs(1:2,:)=[-11,id_fs1,t_departure_fs1,t_arrival_fs1,delv1_fs1,delv2_fs1;
        -12,id_fs2,t_departure_fs2,t_arrival_fs2,delv1_fs2,delv2_fs2];
    star_ID(1,10:11) = [id_fs1 id_fs2]; % Solutions available from J_store
    fprintf(fileID,['\n' repmat('%0.0f,',1,2) repmat('%0.12f,',1,8)],settlement_tree_fs');
elseif id_fs1 ~=0 & id_fs2 ==0
    settlement_tree_fs(1,:)=[-11,id_fs1,t_departure_fs1,t_arrival_fs1,delv1_fs1,delv2_fs1];
    star_ID(1,10) = id_fs1; % Solutions available from J_store
    fprintf(fileID,['\n' repmat('%0.0f,',1,2) repmat('%0.12f,',1,8)],settlement_tree_fs');
elseif id_fs1 ==0 & id_fs2 ~=0
    settlement_tree_fs(1,:)=[-11,id_fs2,t_departure_fs2,t_arrival_fs2,delv1_fs2,delv2_fs2];
    star_ID(1,10) = id_fs2; % Solutions available from J_store
    fprintf(fileID,['\n' repmat('%0.0f,',1,2) repmat('%0.12f,',1,8)],settlement_tree_fs');
end

%% Greedy strategy for settlement
% For departure time of 8 (6+2) myr, identify the three closest stars
% positionally and make a Gen 1 set of settler ships

ss_solver_option=1; % Solver options: 1 for two impulse and 2 for three impulse strategy
gen = 1; % 1st generation of settler ships to be seen
settlement_tree_ss=[];
max_gen=20;

while(gen<max_gen)
    
    tic
    % remove the stars that are already occupied for generating the search
    % space
    disp(['Current evaluation gen:' num2str(gen)])
    
    star_id_settled = star_ID(star_ID~=0); % IDs for settled stars
    
    % Set of upto 3 ^ k stars at Gen k, find 3 closest star for each settled star of Gen k such that upto 3 ^ (k+1) stars selected.
    stars_gen_k = star_ID(gen,:);
    stars_gen_k = stars_gen_k(stars_gen_k~=0);  % find the latest generation settled stars
    
    % Search candidate star pairs for all stars in the gen k
    for  j =1 : length(stars_gen_k)
        
        tof = 2.5;              % starting with small guess
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
            
            indx_current_star = find(settlement_tree_ss(:,2) == ID_jk);     % Find the starID in the arrival ID column of the SS tree
            
            if settlement_tree_ss(indx_current_star(1),3) == 2
                t_departure=settlement_tree_ss(indx_current_star(1),5)+2;       % 2 yr for settling
            else
                t_departure=settlement_tree_ss(indx_current_star(1),6)+2;       % 2 yr for settling
            end
            
            if mod(t_departure,0.5)~=0                                      % In case an arrival time has been calculated by the fmincon solver as a non-multiple of 0.5
                t_departure=t_departure-mod(t_departure,0.5)+0.5;
            end
            
            t_arrival= t_departure+tof;               % myr
            i_arrival = t_arrival/0.5 +1;             % position history column corresponding to the arrival time
            i_departure = t_departure/0.5+1;          % position history column corresponding to the departure time
            
            if t_arrival>89.5
                continue
            end
            
            star_positions_target=[x(2:end,i_arrival),y(2:end,i_arrival),z(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof
        end
        
        
        idx_jk = find(ID_jk == star_data(:,1));                                        % find the row number for the ID in the star database
        r0_jk = [x(idx_jk,i_departure) y(idx_jk,i_departure) z(idx_jk,i_departure)];   % Position of the star ID of the jth element of the kth gen at departure time
        v0_jk = [vx(idx_jk,i_departure) vy(idx_jk,i_departure) vz(idx_jk,i_departure)];% Position of the star ID of the jth element of the kth gen at departure time
        x0 = [r0_jk';v0_jk'];
        
        for l =1:3  % Find transfers to selected 1/2/3 stars with respect to jth settled star of kth generation
            
            %             tof_fit=polyval(poly_n,norm(r0_jk));% guess tof in myr as a multiple of 0.5 myr
            %             tof_fit=tof_fit-mod(tof_fit,0.5);
            tof_fit=2.5;
            t_arrival=t_departure+tof_fit;                      % tof fit
            t_arrival_temp=t_arrival;                           % Initial arrival time
            t_departure_temp=t_departure;
            i_arrival_temp = t_arrival_temp/0.5 +1;             % position history column corresponding to the arrival time
            i_departure_temp = t_departure_temp/0.5+1;          % position history column corresponding to the departure time
            
            ss_solver_option_loop=ss_solver_option;
            
            if t_arrival_temp > 89
                continue
            end
            
            constraints_met =0;
            num_constraint_violation =0;
            
            while constraints_met == 0                          % Loop to find min tof, dV optimum transfer
                
                %                 star_positions_target=[x(2:end,i_arrival_temp),y(2:end,i_arrival_temp),z(2:end,i_arrival_temp)]; % Except sun, all position values for stars at t=tof
                idx_vec = find_closest_momentum_star_ss_strategy6(x0,1,star_ID,x,y,z,i_vec,min_sep,r_max,r_min,min_search_angle,star_database(star_database~=0),t_arrival_temp); % 3 closest stars settlerships
                
                % Update grid center star states
                %                 idx_vec = find_closest_grid_star_strategy6(x0,1,star_ID,x,y,z,t_arrival_temp,star_database(star_database~=0),R_vec); % 3 closest stars settlerships
                
                if isempty(idx_vec)                             % No good unoccupied stars seen
                    break
                end
                
                idx = idx_vec;
                x_t = [x(idx+1,i_arrival_temp) y(idx+1,i_arrival_temp) z(idx+1,i_arrival_temp) vx(idx+1,i_arrival_temp) vy(idx+1,i_arrival_temp) vz(idx+1,i_arrival_temp)]'; % Target states
                
                tof = t_arrival_temp-t_departure_temp;
                r0 = r0_jk';
                rt = x_t(1:3);
                vt = x_t(4:6);
                v0_guess = [vx(idx_jk,i_departure_temp) vy(idx_jk,i_departure_temp) vz(idx_jk,i_departure_temp)]'; % guess value provided as the star velocity itself
                
                
                if ss_solver_option_loop == 1
                    [states,v0,vf]= shooting_solver(r0,rt,tof,v0_guess);
                    dv1= v0-v0_guess'; % transfer impulse
                    delv_transfer = norm(dv1);
                    dv2= -vf+vt';      % rendezvous impulse
                    delv_rendezvous = norm(dv2);
                    
                    dv_consumption = (delv_transfer+delv_rendezvous)* kpcpmyr2kms;  % km/s
                    
                    if norm(dv1) * kpcpmyr2kms > 175 || norm(dv2) * kpcpmyr2kms > 175
                        failure_condition=1;
                    else
                        
                        x0_star_guess=[r0;v0']; % Guess a good solution
                        [t_burn,deltav_ss,exit_flag] = optimal_shooter_hop_ss_strategy5(x0,x_t,tof,x0_star_guess); % contains min time solution
                        
                        t_burn_ss=t_burn+t_departure_temp;
                        
                        cond1_ss= norm(deltav_ss(1,1:3)) * kpcpmyr2kms > 175;
                        cond2_ss= norm(deltav_ss(1,4:6)) * kpcpmyr2kms > 175;
                        cond3_ss= norm(deltav_ss(1,7:9)) * kpcpmyr2kms > 175;
                        cond4_ss= (norm(deltav_ss(1,1:3))+norm(deltav_ss(1,4:6))+norm(deltav_ss(1,7:9)) )* kpcpmyr2kms >400;
                        cond5_ss= exit_flag~=2 && exit_flag~=1;
                        cond6_ss= t_burn_ss(3)>89.5;
                        
                        dv_consumption = (norm(deltav_ss(1,1:3))+norm(deltav_ss(1,4:6))+norm(deltav_ss(1,7:9)) )* kpcpmyr2kms;
                        
                        if  cond1_ss | cond2_ss | cond3_ss | cond4_ss |  cond5_ss | cond6_ss
                            failure_condition=0;    % accept shooter solution
                            dv_consumption = (delv_transfer+delv_rendezvous)* kpcpmyr2kms;  % km/s
                        else 
                            failure_condition=0;
                            ss_solver_option_loop=2;
                        end
                        
                    end
                    
                    
                else
                    disp('Wrong settler ship solver option')
                end
                
                
                if failure_condition==1
                    
                    %Unsuccessful transfer
                    if ss_solver_option_loop==1
                        tof =tof+1;                                         % Increase tof by 0.5
                        t_arrival_temp= t_arrival_temp+1;                   % myr
                        i_arrival_temp=t_arrival_temp/0.5+1;                % column index
                        num_constraint_violation = num_constraint_violation+1;
                    
                    else
                        disp('Wrong option selected for the ss solver')
                    end
                    
                    if num_constraint_violation> 40 || t_arrival_temp>89.5
                        break
                    end
                    
                else
                    
                    % Successful transfer
                    delV_used = delV_used + dv_consumption; % km/s
                    delV_max = delV_max + 400; % km/s
                    error_J_term = J_N_r_theta(star_ID(star_ID~=0),star_data);
                    J_merit = error_J_term * delV_max / delV_used;
                    
                    star_ID(gen+1,(j-1)*3+l)=star_data(idx+1,1);
                    
                    % Settlement Pod line
                    % Write the actions in the txt solution file.
                    
                    if ss_solver_option_loop == 1
                        settlement_tree_ss=[settlement_tree_ss; ID_jk idx 2 t_departure_temp t_arrival_temp dv1* kpcpmyr2kms dv2*kpcpmyr2kms zeros(1,4)];
                        fprintf(fileID,['\n' repmat('%0.0f,', 1, 3) repmat('%0.12f,', 1, 8) ],settlement_tree_ss(end,1:11));
                    elseif ss_solver_option_loop == 2
                        settlement_tree_ss=[settlement_tree_ss; ID_jk idx 3 t_burn_ss deltav_ss*kpcpmyr2kms];
                        fprintf(fileID,['\n' repmat('%0.0f,', 1, 3) repmat('%0.12f,', 1, 12) ],settlement_tree_ss(end,1:15));
                    end
                    
                    ss_solver_option_loop = 1;
                    break
                    
                end
                
            end
            
        end
        
    end
    % Update the departure and arrival times, generation gen
    gen = gen+1;     % next generation of settler ships to be seen
    toc
    toc-tic
end

if J_merit>200
    disp(['J_merit:' num2str(J_merit)])
end


% plot star positions versus generation
for i = 1:13
    for j=1:length(star_ID(star_ID~=0))
        i_f=90/0.5+1;
        if star_ID(i,j)==0
        else
            idx_ij = star_ID(i,j)+1;
            x_t = [x(idx_ij,i_f) y(idx_ij,i_f) z(idx_ij,i_f) vx(idx_ij,i_f) vy(idx_ij,i_f) vz(idx_ij,i_f)]'; % Target states
            figure(1)
            hold on
            plot3(x_t(1),x_t(2),x_t(3),'*')
            %         text(x_t(1),x_t(2),x_t(3),num2str(idx_ij))
            hold on; plot3(0,0,0,'+')
            axis square
            xlim([-32 32])
            ylim([-32 32])
        end
    end
    pause(5)
end
savefig(figname);

