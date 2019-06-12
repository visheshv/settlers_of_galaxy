function [t_burn,deltav,exit_flag] = optimal_shooter_hop_ss_strategy5(star_x0,xt,tof,xi)

% Input 
% star_x0 = initial state of the star as a colmn vector
% xt = final state of the target star as a column vector
% tof = total time of flight of the transfer
% xi =  Initial state of the settler ship ( including modified velocity) 

% Output 
% t_burn: contains t1,t2,t3 (impulse times) myr
% dv: contains dv1x dv1y dv2z dv2x dv2y dv2z dv3x dv3y dv3z (kpc/myr)
% ef: exitflag 2 is good, otherwise bad 

% star_x0 = [6.89029819892678;	-5.96896371360469;	0.0558432694044516; -0.167899193455042;	-0.193813276829153;	0.000198657002247061];
% tof = 19;
% xt = [0.325086955798508; -9.28756499339524; 1.99155694398751; -0.243539204775551; 0.00664123591508843; 0.0707247308801087];

A = []; B = []; Aeq = []; Beq = []; 
options = optimoptions('fmincon','algorithm','sqp','Display','None','MaxFunctionEvaluations',100000,'MaxIterations',60);

%% Now the transfer has to be accomplished with 3 deltaVs and total constant time, broken into 3 steps

tspan = linspace(0,tof,201);
opts = odeset('AbsTol',1e-5);
% x0 = star_x0;
x0=xi;

[integrationtime,states]=ode45(@state_dynamics,tspan,x0,opts);   % containing integrated states

designvariables = vertcat(x0,states(101,:)',tof/2,tof/2);        % containing 14 variables : initial states, intermediate states at tof/2,tof for both segments 

lb = -Inf*ones(length(designvariables),1); lb(13)=1;lb(14)=1;    % set bounds for design variables
ub = Inf*ones(length(designvariables),1);

% run fmincon solver
[X_multiple_hop,cost_multiple_hop,exit_flag,output_multiple_hop] = fmincon(@(x) obj_fn_multiple_hop(x,star_x0,xt,tof), designvariables, A,B,Aeq,Beq,lb,ub,@(x) nonlcon_multiple_hop(x,star_x0,xt,tof),options);


tspan1 = linspace(0,X_multiple_hop(13),100);
opts = odeset('AbsTol',1e-5);
[~,states1]=ode45(@state_dynamics,tspan1,X_multiple_hop(1:6,1),opts);

tspan2 = linspace(0,X_multiple_hop(14),100);
[~,states2]=ode45(@state_dynamics,tspan2,X_multiple_hop(7:12,1),opts);

% Calculate dV vectors for departure, intermediate and arrival deltaVs
dep = ((X_multiple_hop(4:6) - star_x0(4:6)));
inter =  ((-states1(end,4:6)'+X_multiple_hop(10:12)));
arr = ((-states2(end,4:6)' + xt(4:6)));

t_burn = [0,X_multiple_hop(13),X_multiple_hop(14)+X_multiple_hop(13)];
deltav = [dep',inter',arr'];  % dv1, dv2, dv3

end

