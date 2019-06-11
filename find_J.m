idx = find_J(idx,star_ID)

% inputs
% star_ID: existing settled star database
% idx: the indices of the settled stars

star_id_settled = star_ID(star_ID~=0); % IDs for settled stars
star_id_settled = [star_id_settled idx];

R_vec = data(:,2);       % kpc (kiloparsecs)
theta_f_vec = data(:,6); % final polar angle (deg)

n_R = 30;      % number of segments of R_k
n_Theta = 32;  % number of segments of Theta_k

R_k = 2 + (0:n_R);                        % Search radius range kpc
Theta_k = -pi +  (2*pi/32) * (0:n_Theta); % Search angle range radians

s_r = 1;           % kpc
s_theta = 2*pi/32; % radian

x_vec = R_vec * pi/180;
err_r = error_eval(R_k,x_vec,s_r,1);

x_vec = theta_f_vec * pi/180;
err_theta = error_eval(Theta_k,x_vec,s_theta,2);

err_store(count) = err_r+err_theta;
error_J_term(count) = i / (1 + 1e-4 * i * err_store(count));

end