<<<<<<< HEAD
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

num_stars=1e5;
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

theta_f_compute = atan2d(y(:,1),x(:,1));

figure(1)
r=(x.^2+y.^2+z.^2).^0.5;
plot3(x(2:end),y(2:end),z(2:end),'o');
hold on
plot3(x(1),y(1),z(1),'*');

% label(h,strcat('ID: ',num2str(i)),'location','bottom','slope')





=======
% Format of data   
% ID  ,  R (kpc) ,  i (deg)  , Omega (deg),  phi (deg) ,     theta_f (deg)    

clear all
close all

load star_data

kpc2km = 30856775814671900;

k = [0.00287729 0.0023821 -0.0010625 0.000198502 -1.88428e-05 9.70521e-07 -2.70559e-08 3.7516e-10 -1.94316e-12];
R_vec = data(:,2); % kpc (kiloparsecs)
phi_vec = data(:,5); % phi (deg)
omega_vec = data(:,4); % omega (deg)
i_vec = data(:,3); % i (deg)

v_vec = ([R_vec.^0 R_vec.^1 R_vec.^2 R_vec.^3 R_vec.^4 R_vec.^5 R_vec.^6 R_vec.^7 R_vec.^8]*k').^-1; % Star speeds km/s
n_vec = (1/kpc2km)*(v_vec ./ R_vec); % (rad/sec)

myr = 1e6*365*86400;
num_stars=1e5;
span = length(0:myr:90*myr);

x=zeros(num_stars,span); y=zeros(num_stars,span); z=zeros(num_stars,span);
vx=zeros(num_stars,span); vy=zeros(num_stars,span); vz=zeros(num_stars,span);

for i=1:size(v_vec,1)
    counter=0;
    
    for t=0;%:myr:90*myr
        counter=counter+1;
        x(i,counter)= R_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*cosd(omega_vec(i)) - sind(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*sind(omega_vec(i)));
        y(i,counter)= R_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*sind(omega_vec(i)) + sind(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*cosd(omega_vec(i)));
        z(i,counter)= R_vec(i)*(sind(n_vec(i)*t+phi_vec(i))*sind(i_vec(i)));
        vx(i,counter)= v_vec(i)*(-sind(n_vec(i)*t+phi_vec(i))*cosd(omega_vec(i)) - cosd(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*sind(omega_vec(i)));
        vy(i,counter)= v_vec(i)*(-sind(n_vec(i)*t+phi_vec(i))*sind(omega_vec(i)) + cosd(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*cosd(omega_vec(i)));
        vz(i,counter)= v_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*sind(i_vec(i)));
    end
        
end

figure(1)
r=(x.^2+y.^2+z.^2).^0.5;
plot3(x(2:end),y(2:end),z(2:end),'o');
hold on
plot3(x(1),y(1),z(1),'*');

% label(h,strcat('ID: ',num2str(i)),'location','bottom','slope')





>>>>>>> 0aa0ed9179f0d3a0b4129570b10f0c591aa29d32
