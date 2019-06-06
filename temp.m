% Format of data
% ID  ,  R (kpc) ,  i (deg)  , Omega (deg),  phi (deg) ,     theta_f (deg)

clear all
close all

load star_data

kpc2km = 30856775814671900;
myr = 1e6*31557600;
kms2kpcpmyr = myr/kpc2km;
dtr = pi/180;

k = [0.00287729 0.0023821 -0.0010625 0.000198502 -1.88428e-05 9.70521e-07 -2.70559e-08 3.7516e-10 -1.94316e-12];
R_vec = data(:,2); % kpc (kiloparsecs)
phi_vec = data(:,5); % phi (deg)
omega_vec = data(:,4); % omega (deg)
i_vec = data(:,3); % i (deg)
theta_f= data(:,6); % final polar angle (deg)

v_vec = kms2kpcpmyr*([R_vec.^0 R_vec.^1 R_vec.^2 R_vec.^3 R_vec.^4 R_vec.^5 R_vec.^6 R_vec.^7 R_vec.^8]*k').^-1; % Star speeds kpc/Myr
n_vec = (1/dtr)*(v_vec ./ R_vec); % (deg/mYR)

num_stars=1e5+1;
span = length(0:.5:90);

x=zeros(num_stars,span); y=zeros(num_stars,span); z=zeros(num_stars,span);
vx=zeros(num_stars,span); vy=zeros(num_stars,span); vz=zeros(num_stars,span);

for i=1:size(v_vec,1)
    counter=0;
    
    for t=0:0.5:90
        counter=counter+1;
        x(i,counter)= R_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*cosd(omega_vec(i)) - sind(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*sind(omega_vec(i)));
        y(i,counter)= R_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*sind(omega_vec(i)) + sind(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*cosd(omega_vec(i)));
        z(i,counter)= R_vec(i)*(sind(n_vec(i)*t+phi_vec(i))*sind(i_vec(i)));
        vx(i,counter)= v_vec(i)*(-sind(n_vec(i)*t+phi_vec(i))*cosd(omega_vec(i)) - cosd(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*sind(omega_vec(i)));
        vy(i,counter)= v_vec(i)*(-sind(n_vec(i)*t+phi_vec(i))*sind(omega_vec(i)) + cosd(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*cosd(omega_vec(i)));
        vz(i,counter)= v_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*sind(i_vec(i)));
        
    end
       
end


for t=0:0.5:90
    counter=t/0.5+1;
    figure(1)
    plot(x(:,counter),y(:,counter),'.')
    pause(0)
    
end
    

i_temp = find((R_vec(:,1)< (R_vec(1)+0.01)) & (R_vec(:,1)> (R_vec(1)-0.01)) & z );
x0 = [x(1);y(1);z(1);vx(1);vy(1);vz(1)];  % Sol position at t0
x0(4:6) = (norm(x0(4:6))+20 *kms2kpcpmyr)*(x0(4:6)/norm(x0(4:6))); % 30 km/s dV at the start

% Guess input provided [r0;v0]   kpc, kpc/mny
stm0 = zeros(1,36); stm0(1,1)=1; stm0(1,8)=1; stm0(1,15)=1; stm0(1,22)=1; stm0(1,29)=1; stm0(1,36)=1;
tspan = 0:0.5:50;
ic = horzcat(x0',stm0);
[t,states]=ode45(@dynamics,tspan,ic);
% i_temp=i_temp(2);

for k=0:0.5:50
    figure(1)
    r_ms = states(k/0.5+1,1:3);
    v_ms = states(k/0.5+1,4:6);
    %     plot3(r_ms(1),r_ms(2),r_ms(3),'*')
    
    for l=1:length(i_temp)
        %         plot3(x(i_temp(l),k/0.5+1),y(i_temp(l),k/0.5+1),z(i_temp(l),k/0.5+1),'o')
        norm([vx(i_temp(l),k/0.5+1),vy(i_temp(l),k/0.5+1),vz(i_temp(l),k/0.5+1)])
        figure(1)
        hold on
        plot3(vx(i_temp(l),k/0.5+1)-v_ms(1),vy(i_temp(l),k/0.5+1)-v_ms(2),vz(i_temp(l),k/0.5+1)-v_ms(3),'o')
        figure(2)
        hold on
        plot3(x(i_temp(l),k/0.5+1)-r_ms(1),y(i_temp(l),k/0.5+1)-r_ms(2),z(i_temp(l),k/0.5+1)-r_ms(3),'*')
    end
    view(0,90)
%     pause(0.2)
%     close
end


theta_f_compute = atan2d(y(:,end),x(:,end));

% Compare with polar angles given in the data file
err_theta= max(abs(theta_f_compute-theta_f));  % Maximum error value of ~ 1e-4

% save star_snapshots x y z vx vy vz

[R_new,I]=sort(R_vec);
n_vec_new=n_vec(I);
tof_new=(180/16)*(n_vec_new.^-1);
figure(3)
plot(R_new,n_vec_new); hold on;
plot(R_new,tof_new);
xlabel('R kpc');legend('n (deg/myr)','tof (myr)')

%% find the target stars as the center of each grid points at t=0;
r_stars = [x(:,1) y(:,1) z(:,1)];
theta_stars = atan2(r_stars(:,2),r_stars(:,1));
r_stars_norm = vecnorm(r_stars')';
star_target_grids = zeros(30,32);

for radius = 2:1:31 % Move up by 1000 stars to be settled

    for theta = -pi: pi/16:(pi-pi/16)
        
        cond1=r_stars_norm>=radius;
        cond2=r_stars_norm<(radius+1);
        cond3=theta_stars>=theta;
        cond4=theta_stars<(theta+pi/16);
        
        id = find(cond1 & cond2 & cond3 & cond4);
        % find star closest to grid center
        [val,id_center] = min((r_stars_norm(id,:)- (radius+0.5)).^2 + (theta_stars(id,:)- (theta+pi/32)).^2 + r_stars(id,3).^2);
        
        id = id(id_center);
        star_target_grids(radius-1, round((theta+pi)*(16/pi)+1))=id;
        
        hold on
        plot3(x(id,1),y(id,1),z(id,1),'*')
        
    end
end
    


% figure(1)
% r=(x.^2+y.^2+z.^2).^0.5;
% plot3(x(2:end),y(2:end),z(2:end),'o');
% hold on
% plot3(x(1),y(1),z(1),'*');

% label(h,strcat('ID: ',num2str(i)),'location','bottom','slope')




