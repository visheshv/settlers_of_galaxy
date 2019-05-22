%Analysis of shooting method
% Shooting problem

clc
clear;
close all

load star_snapshots
kpc2km = 30856775814671900;
myr = 1e6*31557600;

% Outputs
% analysis of dV vs tof. starting from Sol, find and transfer to the closest
% star at each snap
store_results =zeros(180,19,180);

tic
for j = 1:1:180
disp(['Starting time from Sol : ',num2str((j-1) * 0.5)]);

% Initialize sol state
x0 = [x(1,j) y(1,j) z(1,j) vx(1,j) vy(1,j) vz(1,j)]';  % Sol position at t0
r0 = x0(1:3);
temp_store_results =zeros(180,19);

for i=j+1:181    
    
tof= (i-1) * 0.5;  % Myr

% Sol-less data
star_positions_target=[x(2:100001,i),y(2:100001,i),z(2:100001,i)]; % Except sun, all position values for stars at t=tof
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
tspan = (j-1)*0.5:0.1:tof;
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

% target position, intial velocity, final velocity, target velocity and tof
%store_results(i-1,:) = [x_t(1:6)' ic(1,4:6) states(end, 4:6) x_t(4:6)'  tof]; 
temp_store_results(i-1,:) = [x_t(1:6)' ic(1,1:6) states(end, 1:6) tof]; 
end

store_results(:,:,j) = temp_store_results;
end

toc;
toc-tic

save shooting_analysis_variable_Time

%% Results
%{
start_id = 2;
x0 = [x(1,start_id) y(1,start_id) z(1,start_id) vx(1,start_id) vy(1,start_id) vz(1,start_id)]';

tf=[0:0.5:90]';
delr= vecnorm((store_results(:,1:3,start_id)-x0(1:3,1,1)')')';
delv_transfer= vecnorm((store_results(:,10:12,start_id)-x0(4:6,1,1)')')'*kpc2km/myr;
delv_rendezvous= vecnorm((store_results(:,4:6,start_id)-store_results(:,16:18,start_id))')'*kpc2km/myr;
delv_r_store = [delr delv_transfer delv_rendezvous store_results(:,end,start_id)];

figure(1)
plot(delr,delv_transfer,'o')
hold on
plot(delr,delv_rendezvous,'*')
legend('delv_{transfer} (km/s)','delv_{rendezvous} (km/s)')
xlabel(' |rf-r0| in kpc'); ylabel('delv [km/s]')
grid on;

figure(2)
plot(delv_r_store(:,end),delv_transfer)
ylabel('delv_{transfer} (km/s)')
xlabel('tof (myr)')
grid on;

figure(3)
plot(delv_r_store(:,end), delr)
xlabel('tof (myr)');ylabel('Position of the closest star (kpc)')

%{
figure(4)
plot(tf,delv_rendezvous)
xlabel('tof (myr)');ylabel('dv rendezvous (km/s)')
%}
%}
%% histogram
min_data = zeros(180,4);
for j = 1:1:180
x0 = [x(1,j) y(1,j) z(1,j) vx(1,j) vy(1,j) vz(1,j)]';
sol_time_transfer = (j-1) * 0.5;  
delr= vecnorm((store_results(j:end,1:3,j)-x0(1:3,1,1)')')';
delv_transfer= vecnorm((store_results(j:end,10:12,j)-x0(4:6,1,1)')')'*kpc2km/myr;
delv_rendezvous= vecnorm((store_results(j:end,4:6,j)-store_results(j:end,16:18,j))')'*kpc2km/myr;
%delv_r_store = [delr delv_transfer delv_rendezvous store_results(:,end,j)];

min_data(j,:) = [sol_time_transfer, min(abs(delr)), min(abs(delv_transfer)), min(abs(delv_rendezvous))];

end