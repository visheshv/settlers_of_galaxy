function error_J_term = J_N_r_theta(star_id,data)

% Inputs 
% StarID matrix contains all the starIDs in a vector
% data: star data originally provided

% Outputs: J_term : N / (1 + 1e-4 * N * (E_r + E_theta))

data_settled = data(star_id+1,:);
N = length(data_settled);    

% Lets look at the R_, theta_f data
R_vec = data_settled(:,2);       % kpc (kiloparsecs)
theta_f_vec = data_settled(:,6); % final polar angle (deg)

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

err = err_r+err_theta;
error_J_term = N / (1 + 1e-4 * N * err);
end