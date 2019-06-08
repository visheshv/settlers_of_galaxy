
% Find out possible solutions for star interceptions after the mothership
% has exhausted its last impulse
load star_snapshots.mat
kpc2km = 30856775814671900;
myr = 1e6*31557600;
kms2kpcpmyr = myr/kpc2km;
dtr = pi/180;
options = odeset('Events',@statesEvents);

N = 5000;

fpa1_range = (180*(-1 + 2*rand(1,N)));
fpa2_range = (180*(-1 + 2*rand(1,N)));
fpa3_range = (180*(-1 + 2*rand(1,N)));

del_v1_range = randi(200,N,1) * kms2kpcpmyr;
del_v2_range = randi(200,N,1) * kms2kpcpmyr;
del_v3_range = randi(100,N,1) * kms2kpcpmyr;

t1_range = 1:0.5:20; 
t1_range_indx = randi(length(t1_range),N,1);

t2_range = 1:0.5:10; 
t2_range_indx = randi(length(t2_range),N,1);

t3_range = 1:0.5:10; 
t3_range_indx = randi(length(t3_range),N,1);

t_depart_range = 0:0.5:9.5;
t_depart_indx = randi(length(t_depart_range),N,1);

clear states_log init_cond

states_log = zeros(184,8,N);
init_cond = zeros(N,20);

for ii = 1:1:N
    
    disp(['MC id: ',num2str(ii)]);    
    
    % spiral config: fpa 80,80,120, delv: 200,200,100, t: 20,10,10
    fpa1=fpa1_range(ii);
    fpa2=fpa2_range(ii);
    fpa3=fpa3_range(ii);
    
    del_v1 = del_v1_range(ii);
    del_v2 = del_v2_range(ii);
    del_v3 = del_v3_range(ii);
    
    t1=t1_range(t1_range_indx(ii));
    t2=t2_range(t2_range_indx(ii));
    t3=t3_range(t3_range_indx(ii));
    
    t_depart = t_depart_range(t_depart_indx(ii));
    i_query = t_depart/0.5 + 1;
    
    r_mat1=[cosd(fpa1) sind(fpa1) 0; -sind(fpa1) cosd(fpa1) 0; 0 0 1];
    r_mat2=[cosd(fpa2) sind(fpa2) 0; -sind(fpa2) cosd(fpa2) 0; 0 0 1];
    r_mat3=[cosd(fpa3) sind(fpa3) 0; -sind(fpa3) cosd(fpa3) 0; 0 0 1];
    
    r0=[x(1,i_query) y(1,i_query) z(1,i_query)]';
    v0=[vx(1,i_query) vy(1,i_query) vz(1,i_query)]';
    v0=v0+del_v1 *r_mat1 *v0/norm(v0);
    
    x0 = [r0;v0];  % Sol position at t0
    stm0 = zeros(1,36); stm0(1,1)=1; stm0(1,8)=1; stm0(1,15)=1; stm0(1,22)=1; stm0(1,29)=1; stm0(1,36)=1;
    tspan = 0:0.5:t1;
    ic = horzcat(x0',stm0);
    [time,states,TE,VE]=ode45(@dynamics,tspan,ic,options);
    
    r0=states(end,1:3)';v0=states(end,4:6)'+ del_v2 *r_mat2 * states(end,4:6)'/norm(states(end,4:6)');x0 = [r0;v0];  % Sol position at t0
    tspan = 0:0.5:t2;
    ic = horzcat(x0',stm0);
    [time1,states1,TE1,VE1]=ode45(@dynamics,tspan,ic,options);
    
    r0=states1(end,1:3)';v0=states1(end,4:6)' + del_v3 *r_mat3* states1(end,4:6)'/norm(states1(end,4:6)'); x0 = [r0;v0];  % Sol position at t0
    tspan = 0:0.5:t3;
    ic = horzcat(x0',stm0);
    [time2,states2,TE2,VE2]=ode45(@dynamics,tspan,ic,options);
    
    r0=states2(end,1:3)';v0=states2(end,4:6)'; x0 = [r0;v0];  % Sol position at t0
    tspan = 0:0.5:90-(t1+t2+t3);
    ic = horzcat(x0',stm0);
    [time3,states3,TE3,VE3]=ode45(@dynamics,tspan,ic,options);
    
    r_hist=[states(:,1:3);states1(:,1:3);states2(:,1:3);states3(:,1:3)];
    v_hist=[states(:,4:6);states1(:,4:6);states2(:,4:6);states3(:,4:6)];
    t_hist = [time; t1+time1; (t1+t2) + time2; (t1+t2+t3) + time3];
    event_log = zeros(length(t_hist),1);
    event_count = 0;
    
    if ~isempty(TE)
        for jj = 1:1:length(TE)
            idx = find(abs(t_hist-TE(jj))<10e-3);
            event_log(idx(1),1) = 1;
            event_count = event_count + 1;
        end
    else
        TE = inf;
    end
    
    if ~isempty(TE1)
        for jj = 1:1:length(TE1)
            idx1 = find(abs(t_hist-(TE1(jj) + t1))<10e-3);
            event_log(idx1(1),1) = 1;
            event_count = event_count + 1;
        end
    else
        TE1 = inf;
    end
    
    if ~isempty(TE2)
        for jj = 1:1:length(TE2)
            idx2 = find(abs(t_hist-(TE2(jj)+t1+t2))<10e-3);
            event_log(idx2(1),1) = 1;
            event_count = event_count + 1;
        end
    else
        TE2 = inf;
    end
    
    if ~isempty(TE3)
        for jj = 1:1:length(TE3)
            idx3 = find(abs(t_hist-(TE3(jj)+t1+t2+t3))<10e-3);
            event_log(idx3(1),1) = 1;
            event_count = event_count + 1;
        end
    else
        TE3 = inf;
    end
    
    init_data = [fpa1 fpa2 fpa3 [del_v1 del_v2 del_v3]*(1/kms2kpcpmyr) t1 t2 t3 i_query event_count TE t1+TE1' (t1+t2+TE2') (t1+t2+t3+TE3')];
    init_cond(ii,:) = [init_data , zeros(1,20-length(init_data))];
    states_log(:,:,ii) = [t_hist r_hist v_hist event_log];
    
end
 
save mothership_mc.mat init_cond states_log;

%% failure:
success_idx = find(init_cond(:,11) == 0);
failure_idx = find(init_cond(:,11) ~= 0);

for kk = 1:1:length(failure_idx)
plot3(states_log(:,2,failure_idx(kk)),states_log(:,3,failure_idx(kk)),states_log(:,4,failure_idx(kk))); hold on;
end


hist(init_cond(success_idx,1));


%% plot
for kk = 1:1:N
plot3(states_log(:,2,kk),states_log(:,3,kk),states_log(:,4,kk)); hold on;
end