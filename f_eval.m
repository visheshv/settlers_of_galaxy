function f_val = f_eval(X_k,x_vec,s)

% Evaluate the f function 
% Input
% X_k: Query state, e.g. R_k,Theta_k 
% x_vec: Vector of the states of settled stars, e.g. r_i or theta_i
% s: Threshold state position or angle
%
% Output
% f function

N=length(x_vec); % number of settled stars
f_val=0;         % Define intial value

for i = 1:N
    
    xi = x_vec(i);
    
    if abs(X_k - xi) >= s
        f = 0;
    else
        f = 1/s - abs(X_k-xi)/s^2; 
    end
    
    f_val = f_val + f;
    
end

f_val = f_val / N;