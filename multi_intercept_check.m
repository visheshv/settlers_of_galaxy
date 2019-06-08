function [store_close_stars]=multi_intercept_check(t_hist,r_hist,v_hist,R_vec,phi_vec,omega_vec,i_vec,v_vec,n_vec)

store_close_stars=zeros(length(t_hist),2);

for j=1:length(t_hist)

    r=r_hist(j,:);
    t=t_hist(j);
    v=v_hist(j,:);
    
    star_positions_target=zeros(1e5+1,3);
    star_vel_target=zeros(1e5+1,3);
    
    for i=1:length(R_vec)
        
    x= R_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*cosd(omega_vec(i)) - sind(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*sind(omega_vec(i)));
    y= R_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*sind(omega_vec(i)) + sind(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*cosd(omega_vec(i)));
    z= R_vec(i)*(sind(n_vec(i)*t+phi_vec(i))*sind(i_vec(i)));
    vx= v_vec(i)*(-sind(n_vec(i)*t+phi_vec(i))*cosd(omega_vec(i)) - cosd(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*sind(omega_vec(i)));
    vy= v_vec(i)*(-sind(n_vec(i)*t+phi_vec(i))*sind(omega_vec(i)) + cosd(n_vec(i)*t+phi_vec(i))*cosd(i_vec(i))*cosd(omega_vec(i)));
    vz= v_vec(i)*(cosd(n_vec(i)*t+phi_vec(i))*sind(i_vec(i)));
    
    if ((i_vec(i) <=180) && (i_vec(i) >=177) || (i_vec(i) >=0) && (i_vec(i) <=3))
        if i==1
            star_positions_target(i,:)=[inf inf inf];
            star_vel_target(i,:)=[inf inf inf];
        else
            star_positions_target(i,:)=[x y z];
            star_vel_target(i,:)=[vx vy vz];
        end
    else
        star_positions_target(i,:)=[inf inf inf];
        star_vel_target(i,:)=[inf inf inf];
    end
        
    end
    
    idx = knnsearch(star_positions_target,r,'K',1);     
    dist_min=norm(star_positions_target(idx,:)- r);
    vel_rendezvous = norm(-v + star_vel_target(idx,:));
    kpc2km = 30856775814671900;
    myr = 1e6*31557600;
    kms2kpcpmyr = myr/kpc2km;
    
    if vel_rendezvous <300 *kms2kpcpmyr
    store_close_stars(j,:)= [idx-1 dist_min];  % -1 to identify the star ID
    end
    
end

    

