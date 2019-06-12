% Shooting solver
% 
% Input 
% r0: Initial position (kpc) 3 x 1
% rt: Target position (kpc) 
% tof: time of flight myr
% v0_guess: Initial guess for velocity at initial time 3 x 1
%
% Output
% [r(t),v(t)]: States
% vi: Initial velocity (kpc/myr)
% vf: Final velocity (kpc/myr)


function [states,v0,vf]= shooting_solver(r0,rt,tof,v0_guess)

x0 = [r0;v0_guess];  % Sol position at t0

% Guess input provided [r0;v0]   kpc, kpc/mny
stm0 = zeros(1,36); stm0(1,1)=1; stm0(1,8)=1; stm0(1,15)=1; stm0(1,22)=1; stm0(1,29)=1; stm0(1,36)=1;

% Propagate dynamics and find X(tf)
tspan = 0:0.1:tof;
del_xf=[1;0;0];
ic = horzcat(x0',stm0);

while_count=0;

while norm(del_xf) > 10^-10
    
    while_count=while_count+1;
    
    % norm(del_xf)
    opts = odeset('AbsTol',1e-8,'RelTol',1e-8);
    [~,states]=ode45(@dynamics,tspan,ic,opts);
    
    del_xf = -[states(end,1)-rt(1);states(end,2)-rt(2);states(end,3)-rt(3)];
    stm_f = [states(end,7:12);states(end,13:18);states(end,19:24);states(end,25:30);states(31:36);states(end,37:42)];
    rel_stm_f = stm_f(1:3,4:6);
    del_v0 = rel_stm_f\del_xf;
    ic(1,4:6)=ic(1,4:6)+del_v0';
    
%     plot3(states(:,1),states(:,2),states(:,3))
%     hold on   
    
if while_count>25
    states(1,4:6)=[inf inf inf];
    states(end,4:6)=[inf inf inf];
    tof=inf;
    break
end

end

% target position, intial velocity, final velocity, target velocity and tof
v0 = states(1,4:6);
vf = states(end,4:6);

end
