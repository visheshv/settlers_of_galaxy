function [store_close_stars]=one_star_targeting_mothership(t_hist,r_hist)

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

store_close_stars=zeros(length(t_hist),2);

for j=1:length(t_hist)
    r=r_hist(j);
    t=t_hist(j);
    
    for i=1:length(R_vec)
        
    x= R_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*cosd(omega_vec(i)) - sind(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*sind(omega_vec(i)));
    y= R_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*sind(omega_vec(i)) + sind(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*cosd(omega_vec(i)));
    z= R_vec(i)*(sind(n_vec(i)*t+phi_vec(i))*sind(i_vec(i)));
    vx= v_vec(i)*(-sind(n_vec(i)*t+phi_vec(i))*cosd(omega_vec(i)) - cosd(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*sind(omega_vec(i)));
    vy= v_vec(i)*(-sind(n_vec(i)*t+phi_vec(i))*sind(omega_vec(i)) + cosd(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*cosd(omega_vec(i)));
    vz= v_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*sind(i_vec(i)));
    
    if i_vec(i) <3
        star_positions_target(i,:)=[x y z];
    else
        star_positions_target(i,:)=[inf inf inf];
    end
        
    end
    
    idx = knnsearch(star_positions_target,r,'K',1);     
    dist_min=norm(star_positions_target(idx,:)- r);
    store_close_stars(j,:)= [idx dist_min];
    
end

    

