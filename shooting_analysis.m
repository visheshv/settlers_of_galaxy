% Analysis of shooting method
% Shooting problem

clc
clear all
close all

load star_snapshots

% Outputs
% analysis of dV vs tof. starting from Sol, find and transfer to the closest
% star at each snap

% Initialize sol state
x0 = [x(1,1) y(1,1) z(1,1) vx(1,1) vy(1,1) vz(1,1)]';  % Sol position at t0
r0 = x0(1:3);

store_results =zeros(180,13);
tic

for i=2:181

tof= (i-1) * 0.5  % Myr

% Sol-less data
star_positions_target=[x(2:end,i),y(2:end,i),z(2:end,i)]; % Except sun, all position values for stars at t=tof
idx = knnsearch(star_positions_target,x0(1:3)'); % Closest star

x_t = [x(idx+1,i) y(idx+1,i) z(idx+1,i) vx(idx+1,i) vy(idx+1,i) vz(idx+1,i)]';

% Guess input provided [r0;v0]   kpc, kpc/mny
stm0 = zeros(1,36);
stm0(1,1)=1;
stm0(1,8)=1;
stm0(1,15)=1;
stm0(1,22)=1;
stm0(1,29)=1;
stm0(1,36)=1;

% Propagate dynamics and find X(tf)
tspan = 0:0.1:tof;
del_xf=[1;0;0];
ic = horzcat(x0',stm0);

while_count=0;

while norm(del_xf) > 10^-10
    
    while_count=while_count+1;
    
    % norm(del_xf)
    [t,states]=ode45(@dynamics,tspan,ic);
    
    del_xf = -[states(end,1)-x_t(1);states(end,2)-x_t(2);states(end,3)-x_t(3)];
    stm_f = [states(end,7:12);states(end,13:18);states(end,19:24);states(end,25:30);states(31:36);states(end,37:42)];
    rel_stm_f = stm_f(1:3,4:6);
    del_v0 = rel_stm_f\del_xf;
    ic(1,4:6)=ic(1,4:6)+del_v0';
    
%     plot3(states(:,1),states(:,2),states(:,3))
%     hold on
if while_count>30
    tof=inf;
end

end

store_results(i-1,:) = [x_t(1:3)' ic(1,4:6) states(end, 4:6) x_t(4:6)'  tof]; % target position, 

end

toc;
toc-tic


tf=0.5:0.5:90;
delr= vecnorm((store_results(:,1:3)-repmat(x0(1:3)',180,1))')';
delv_transfer= vecnorm((store_results(:,4:6)-store_results(:,7:9))')';
delv_rendezvous= vecnorm((store_results(:,7:9)-store_results(:,10:12))')';
delv_r_store = [delr delv_transfer delv_rendezvous store_results(:,end)];

figure(1)
plot(delr,delv_transfer,'o')
hold on
plot(delr,delv_rendezvous,'*')
legend('delv_{transfer} (kpc/myr)','delv_{rendezvous} (kpc/myr)')
xlabel(' |rf-r0| in kpc')


figure(2)
plot(delv_transfer,delv_r_store(:,end))
xlabel('delv_{transfer} (kpc/myr)')
ylabel('tof (myr)')

figure(3)
plot(tf, delr)
xlabel('tof (myr)');ylabel('Position of the closest star (kpc)')

figure(4)
plot(tf,delv_transfer)
xlabel('tof (myr)');ylabel('transfer dV for the closest star (kpc/myr)')



% Propagate STM over the X(t) and find STM(tf), calculate B (3 x 3) matrix





% Update the guess using the B matrix
