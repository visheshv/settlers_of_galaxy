% Find out possible solutions for star interceptions after the mothership
% has exhausted its last impulse
load star_snapshots.mat
kpc2km = 30856775814671900;
myr = 1e6*31557600;
kms2kpcpmyr = myr/kpc2km;
dtr = pi/180;

% spiral config: fpa 80,80,120, delv: 200,200,100, t: 20,10,10

fpa1=80;
fpa2=80;
fpa3=120;

r_mat1=[cosd(fpa1) sind(fpa1) 0; -sind(fpa1) cosd(fpa1) 0; 0 0 1];
r_mat2=[cosd(fpa2) sind(fpa2) 0; -sind(fpa2) cosd(fpa2) 0; 0 0 1];
r_mat3=[cosd(fpa3) sind(fpa3) 0; -sind(fpa3) cosd(fpa3) 0; 0 0 1];

del_v1=200 *kms2kpcpmyr;
del_v2=200 *kms2kpcpmyr;
del_v3=100 *kms2kpcpmyr;

t1=20;
t2=10;
t3=10;

i_query=1;

r0=[x(1,i_query) y(1,i_query) z(1,i_query)]';
v0=[vx(1,i_query) vy(1,i_query) vz(1,i_query)]';
v0=v0+del_v1 *r_mat1 *v0/norm(v0);

x0 = [r0;v0];  % Sol position at t0
stm0 = zeros(1,36); stm0(1,1)=1; stm0(1,8)=1; stm0(1,15)=1; stm0(1,22)=1; stm0(1,29)=1; stm0(1,36)=1;
tspan = 0:0.1:t1;
ic = horzcat(x0',stm0);
[t,states]=ode45(@dynamics,tspan,ic);

r0=states(end,1:3)';v0=states(end,4:6)'+ del_v2 *r_mat2 * states(end,4:6)'/norm(states(end,4:6)');x0 = [r0;v0];  % Sol position at t0
tspan = 0:0.1:t2;
ic = horzcat(x0',stm0);
[t,states1]=ode45(@dynamics,tspan,ic);

r0=states1(end,1:3)';v0=states1(end,4:6)' + del_v3 *r_mat3* states1(end,4:6)'/norm(states1(end,4:6)'); x0 = [r0;v0];  % Sol position at t0
tspan = 0:0.1:t3;
ic = horzcat(x0',stm0);
[t,states2]=ode45(@dynamics,tspan,ic);

r0=states2(end,1:3)';v0=states2(end,4:6)'; x0 = [r0;v0];  % Sol position at t0
tspan = 0:0.1:90-(t1+t2+t3);
ic = horzcat(x0',stm0);
[t,states3]=ode45(@dynamics,tspan,ic);

r_hist=[states(:,1:3);states1(:,1:3);states2(:,1:3);states3(:,1:3)];
v_hist=[states(:,4:6);states1(:,4:6);states2(:,4:6);states3(:,4:6)];

figure(1)
plot3(r_hist(:,1),r_hist(:,2),r_hist(:,3))
xlim([-32 32])
ylim([-32 32])

vel=vecnorm(v_hist')';
figure(2)
plot(vel)
ylabel('vel (kpc/myr)')
