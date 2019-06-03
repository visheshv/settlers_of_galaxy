% Implement a simple strategy to populate stars only with three motherships

% 3. Find the settlement tree

clc
clear all
close all

%% 
% load strat2data
%%{
% Load data
load star_snapshots
load star_data
load J_store
load star_target_grids
 
global star_data  R_vec phi_vec omega_vec i_vec theta_f
global v_vec n_vec
global J_merit delV_max delV_used
global R_min R_max
global x y z vx vy vz

warning('off','all')

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

J_merit =0;
delV_max = 0;
delV_used = 0;
J_store=unique(J_store,'rows');
J_store=sortrows(J_store,6,'descend');
star_ID(1,10:11) = [J_store(1,1) J_store(3,1)]; % Solutions available from J_store
% idx_ms=[9;16;24]; % rows 9,16,24 represent best compromise between max merit value and min time of arrival so that the settling pods have enough time to breed
idx_ms=[283;319;207]; % rows 283,319,207 represent min time of arrival so that the settling pods have enough time to breed
fileID = fopen('strategy3.txt','w');
fprintf(fileID,'%s','strategy3');

settlement_tree_sp=[];
settlement_tree_ms=[];

%% Mother ship capture of near star with minimum transfer deltaV
% Initialize sol state
for mothership_id=-1:-1:-3

    t_departure=J_store(idx_ms(-mothership_id),2);     %myr
    i_departure=t_departure/0.5+1;              %column number query
    t_arrival= J_store(idx_ms(-mothership_id),4);      %myr
    i_arrival=t_arrival/0.5+1;                  %column index corresponding to t = 6myr
    x0 = [x(1,i_departure) y(1,i_departure) z(1,i_departure) vx(1,i_departure) vy(1,i_departure) vz(1,i_departure)]';  % Sol position at t0
    r0 = x0(1:3);
    v0_guess = x0(4:6); % Guess initial condition for solver

    num_impulses_ms=0;
    sum_impulses_ms=0;
    store_results =zeros(1,16);
    delv_store=[];
    state_store=[];

    J_merit_temp =0;
    delV_max_temp =0;
    delV_used_temp =0;

    tic
    constraints_met=0; % check whether delta V constraints are met after finding a solution
    num_violations=0;
    %%

    while constraints_met==0

        % Sol-less data
        if num_impulses_ms==0
            idx = J_store(idx_ms(-mothership_id),1);  % Targets queried from good Sol transfers from J_store
            is_bad_solution =0;
        else            
            [idx, t_departure,t_arrival,is_bad_solution,x0_departure] = find_grid_stars(x0,1,star_ID, t_min_departure,1,inf,J_merit,delV_max,delV_used,x,y,z,vx,vy,vz,star_data,star_target_grids); % Mothership position propagated 1.5 myr since the previous star interception
            r0=x0_departure(1:3);v0_guess= x0_departure(4:6);
        end

        x_t = [x(idx+1,i_arrival) y(idx+1,i_arrival) z(idx+1,i_arrival) vx(idx+1,i_arrival) vy(idx+1,i_arrival) vz(idx+1,i_arrival)]'; % Target states
        rt = x_t(1:3);vt = x_t(4:6);
        
        tof=t_arrival-t_departure;

        [states,v0,vf]= shooting_solver(r0,rt,tof,v0_guess); % Run the algorithm
        store_results = [x_t(1:6)' states(1,4:6) states(end, 4:6) x_t(4:6)'  tof]; % State Store

        delv_transfer= norm(v0-v0_guess');
        delv1= v0-v0_guess';
        delv_rendezvous= vecnorm((store_results(:,10:12)-store_results(:,13:15))')';
        delv2= (-store_results(:,10:12)+store_results(:,13:15)); % rendezvous delv


        %% Identify successful or failed transfer
        if delv_transfer * kpcpmyr2kms > 200 || delv_rendezvous * kpcpmyr2kms > 300 || is_bad_solution ~= 0 % If Mothership initial impulse exceeds  200 km/s or setlling pod impulse exceeds 300 km/s
%             disp(['no solution at tof (myr):' num2str(tof)])
            num_violations=num_violations+1;
            if num_violations>10
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
                % Generation 0 mothership star captures and update the star
                % settled data base
                %             star_database_ms = find_solution_intercepts_ms(t_departure,r0,v0_guess,idx_prev);
                break
            elseif num_impulses_ms==1
                star_ID(1,num_impulses_ms * -mothership_id) = idx; % 1st settlement from Mothership
                delV_max = delV_max + 500; % Add the mothership dv_limits (km/s)
                error_J_term = J_N_r_theta(star_ID(star_ID~=0),star_data);
                J_merit = error_J_term * delV_max / delV_used;
            else
                star_ID(1,num_impulses_ms* -mothership_id) = idx; % nth settlement from Mothership
                error_J_term = J_N_r_theta(star_ID(star_ID~=0),star_data);
                J_merit = error_J_term * delV_max / delV_used;
            end

            delv_store = [delv_store; t_departure delv1* kpcpmyr2kms t_arrival delv2* kpcpmyr2kms];
            settlement_tree_sp =[settlement_tree_sp; mothership_id num_impulses_ms idx t_arrival delv2* kpcpmyr2kms];

            % Update mothership states
            t_min_departure=t_arrival+1.5;                %myr since maneuver has to happen >1 yr after SP rendezvous maneuver

            % Propagate to t_departure
            stm0 = zeros(1,36); stm0(1,1)=1; stm0(1,8)=1; stm0(1,15)=1; stm0(1,22)=1; stm0(1,29)=1; stm0(1,36)=1;
            tspan = 0:0.1:1.5;
            ic = horzcat(states(end,1:6),stm0);
            opts = odeset('AbsTol',1e-5);
            [~,states]=ode45(@dynamics,tspan,ic,opts);

            x0 = states(end,1:6)';                     %state
            r0 =  x0(1:3);
            v0_guess = x0(4:6); % Guess initial condition for solver
            idx_prev = idx;
        end

    end
    %% Mothership line
    length_ms_info=size(delv_store,1);
    dvrow=reshape(delv_store(:,2:4)',1,numel(delv_store(:,2:4)));
    settlement_tree_ms=[mothership_id, 0, length_ms_info, length_ms_info, delv_store(:,1)',dvrow];
    fprintf(fileID,['\n' repmat('%d,',1,4) repmat('%0.12f,',1,4*length_ms_info)],settlement_tree_ms(end,:)');
    % Settlement Pod line
    fprintf(fileID,['\n' repmat('%0.0f,',1, 3) repmat('%0.12f,',1,4)],settlement_tree_sp(end-length_ms_info+1:end,:)');
end

toc;
toc-tic
%% Fastship
% Fastship line

settlement_tree_fs(1:2,:)=[-11,76457,8,61,9.82419017275941,-12.2229794464399,-0.758833716987316,-7.05705080167052,12.9841046519746,-0.931972851585255;
                           -12,46609,2.5,57,-1.51043054136710,-8.21772704567083,5.21514694136725,-0.366263344844732,12.7921582425029,-11.5937056669336];  % row 1 and row3 of initial solution
fprintf(fileID,['\n' repmat('%0.0f,',1,2) repmat('%0.12f,',1,8)],settlement_tree_fs');


save strat2data
%}

%% Greedy strategy for settlement
% For departure time of 8 (6+2) myr, identify the three closest stars
% positionally and make a Gen 1 set of settler ships

gen = 1; % 1st generation of settler ships to be seen
settlement_tree_ss=[];

while(gen<13)
    
    % remove the stars that are already occupied for generating the search
    % space
    disp(['Current evaluation gen:' num2str(gen)])
    
    idx_search = star_data(2:end,1);       % Complete search space for stars except sol
    star_id_settled = star_ID(star_ID~=0); % IDs for settled stars
        
    % Set of upto 3 ^ k stars at Gen k, find 3 closest star for each settled
    % star of Gen k such that upto 3 ^ (k+1) stars selected.
    stars_gen_k = star_ID(gen,:);
    stars_gen_k = stars_gen_k(stars_gen_k~=0); % find the latest generation settled stars
    %r_temp = star_positions_target;
    
    % Search star pairs for all stars in the gen k
    for  j =1 : length(stars_gen_k)
        
        ID_jk = stars_gen_k(j); % star ID of the jth element of the kth gen
        
        if gen==1
            
            if j<=(length(stars_gen_k)-2)
                t_min_departure=settlement_tree_sp(j,4)+2;    % 2 yr for settling a star system
            else
                t_min_departure=settlement_tree_fs(j+2-length(stars_gen_k),4)+2;    % 2 yr for settling a star system
            end
            
        else
            
            indx_current_star = find(settlement_tree_ss(:,2) == ID_jk);
            t_min_departure=settlement_tree_ss(indx_current_star,5)+2; % 2 yr for settling a star system
            
        end
        
        if t_min_departure>=89.5 % Dont execute below this line for the loop if this condition comes
                continue 
        end
                
        x0=[x(ID_jk+1,t_min_departure/0.5+1),y(ID_jk+1,t_min_departure/0.5+1),z(ID_jk+1,t_min_departure/0.5+1),vx(ID_jk+1,t_min_departure/0.5+1),vy(ID_jk+1,t_min_departure/0.5+1),vz(ID_jk+1,t_min_departure/0.5+1)]';
        
        [idx_vec, t_departure,t_arrival,is_bad_solution, ~] = find_grid_stars(x0,3,star_ID, t_min_departure,2,ID_jk,J_merit,delV_max,delV_used,x,y,z,vx,vy,vz,star_data,star_target_grids);
        
        star_positions_target=[x(2:end,i_arrival),y(2:end,i_arrival),z(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof
        star_velocities_target=[vx(2:end,i_arrival),vy(2:end,i_arrival),vz(2:end,i_arrival)]; % Except sun, all position values for stars at t=tof
        
        idx_jk = find(ID_jk == star_data(:,1)); % find the row number for the ID in the star database
        
        if is_bad_solution ~= 3
            
            for l =1:(3-is_bad_solution)  % Find all the closest 3 stars with respect to jth settled star of kth generation
                
                i_departure=t_departure(l)/0.5+1;
                r0_jk = [x(idx_jk,i_departure) y(idx_jk,i_departure) z(idx_jk,i_departure)];   % Position of the star ID of the jth element of the kth gen at departure time
                v0_jk = [vx(idx_jk,i_departure) vy(idx_jk,i_departure) vz(idx_jk,i_departure)];   % Position of the star ID of the jth element of the kth gen at departure time
                x0 = [r0_jk';v0_jk'];
                
                t_arrival_temp=t_arrival(l);
                t_departure_temp=t_departure(l);
                i_arrival_temp = t_arrival_temp/0.5 +1;             % position history column corresponding to the arrival time
                i_departure_temp = t_departure_temp/0.5+1;          % position history column corresponding to the departure time
                
                constraints_met =0;
                num_constraint_violation =0;
                
                idx = idx_vec(l);
                x_t = [x(idx+1,i_arrival_temp) y(idx+1,i_arrival_temp) z(idx+1,i_arrival_temp) vx(idx+1,i_arrival_temp) vy(idx+1,i_arrival_temp) vz(idx+1,i_arrival_temp)]'; % Target states
                
                % Solve the shooting problem for all star pairs and remove the pairs which
                % are problematic in terms of transfer deltaV, rendezvous deltaV,  tof =inf
                tof = t_arrival_temp-t_departure_temp;
                r0 = r0_jk';
                rt = x_t(1:3);
                vt = x_t(4:6);
                v0_guess = v0_jk'; % guess value provided as the star velocity itself
                
                [states,v0,vf]= shooting_solver(r0,rt,tof,v0_guess);
                
                dv1= v0-v0_guess'; % transfer impulse
                delv_transfer = norm(dv1);
                dv2= -vf+vt';      % rendezvous impulse
                delv_rendezvous = norm(dv2);
                
                if norm(dv1) * kpcpmyr2kms > 175 || norm(dv2) * kpcpmyr2kms > 175
                    
                    num_constraint_violation = num_constraint_violation+1;                    
                    if num_constraint_violation> 20 || t_arrival_temp>89.5
                        break
                    end
                    
                else
                    
                    % Successful transfer
                    delV_used = delV_used + delv_transfer * kpcpmyr2kms + delv_rendezvous * kpcpmyr2kms; % km/s
                    delV_max = delV_max + 400; % km/s
                    error_J_term = J_N_r_theta(star_ID(star_ID~=0),star_data);
                    J_merit = error_J_term * delV_max / delV_used;
                    
                    check_star = 1-sum(star_ID(gen+1,:)==idx); % Goes to 1 if idx is a new unique star in the current generation selection
                    
                    if check_star ==1
                        star_ID(gen+1,(j-1)*3+l)=idx;
                        % Settlement Pod line
                        % Write the actions in the txt solution file.
                        settlement_tree_ss=[settlement_tree_ss; ID_jk idx 2 t_departure_temp t_arrival_temp dv1* kpcpmyr2kms dv2* kpcpmyr2kms];
                        fprintf(fileID,['\n' repmat('%0.0f,', 1, 3) repmat('%0.12f,', 1, 8) ],settlement_tree_ss(end,:)');
                    end
                end
            end
        else
        end
    end
    
    % Update the departure and arrival times, generation gen
    gen = gen+1;     % next generation of settler ships to be seen
    
end

% plot star positions versus generation
for i = 1:9
    for j=1:3^(i-1)
        i_f=90/0.5+1;
        if star_ID(i,j)==0
        else
            temp = find(star_ID(i,j) == star_data(:,1)); % find the row number for the ID in the star database
            idx_ij =  star_data(temp,1);
            x_t = [x(idx_ij,i_f) y(idx_ij,i_f) z(idx_ij,i_f) vx(idx_ij,i_f) vy(idx_ij,i_f) vz(idx_ij,i_f)]'; % Target states
            figure(1)
            hold on
            plot3(x_t(1),x_t(2),x_t(3),'*')
            
        end
        
    end
    pause(1)
    
end

