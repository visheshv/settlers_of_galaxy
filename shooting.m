% Shooting problem example

% Inputs
% r0
% rf 
% v0
% tf
clc
clear all
close all

% Outputs
% v0 , vf such there is a trajectory solution between r0,rf
tic
% Simulated inputs
tof= 100;
x_f= [3.1340;3.5442;0.1294;0.1793;-0.1619;0.0906];
x_t = [-6,-10,0];
% Guess input provided [r0;v0]   kpc, kpc/mny
stm0 = zeros(1,36);
stm0(1,1)=1;
stm0(1,8)=1;
stm0(1,15)=1;
stm0(1,22)=1;
stm0(1,29)=1;
stm0(1,36)=1;

x0 = [8.34,0,0,0,-0.2626,0];

% Propagate dynamics and find X(tf)
tspan = 0:0.1:tof;
del_xf=[1;0;0];
ic = horzcat(x0,stm0);

while norm(del_xf) > 10^-10
 norm(del_xf)   
[t,states]=ode45(@dynamics,tspan,ic);

del_xf = -[states(end,1)-x_t(1);states(end,2)-x_t(2);states(end,3)-x_t(3)];
stm_f = [states(end,7:12);states(end,13:18);states(end,19:24);states(end,25:30);states(31:36);states(end,37:42)];
rel_stm_f = stm_f(1:3,4:6);
del_v0 = rel_stm_f\del_xf;
ic(1,4:6)=ic(1,4:6)+del_v0';
plot3(states(:,1),states(:,2),states(:,3))
hold on   
end

toc
% Propagate STM over the X(t) and find STM(tf), calculate B (3 x 3) matrix

plot3(states(:,1),states(:,2),states(:,3))



% Update the guess using the B matrix
