function error_J_term = J_analysis_settled_stars(star_id,data)

data_settled = data(star_id+1,:);

% Sort data by radial position: Closest stars first
data_sorted = sortrows(data_settled,2); % Sort ascending by radial position / theta position 2/6 
N = length(data_settled);    
data = data_sorted(1:N,:); % filter the first 'i' stars

% Lets look at the R_, theta_f data
R_vec = data(:,2);       % kpc (kiloparsecs)
theta_f_vec = data(:,6); % final polar angle (deg)

n_R = 30;      % number of segments of R_k
n_Theta = 32;  % number of segments of Theta_k

R_k = 2 + (0:n_R);                        % Search radius range kpc
Theta_k = -pi +  (2*pi/32) * (0:n_Theta); % Search angle range radians

s_r = 1;           % kpc
s_theta = 2*pi/32; % radian

x_vec = R_vec;
err_r = error_eval(R_k,x_vec,s_r,1);

x_vec = theta_f_vec * pi/180;
err_theta = error_eval(Theta_k,x_vec,s_theta,2);

err_store = err_r+err_theta;
error_J_term = N / (1 + 1e-4 * N * err_store);
end