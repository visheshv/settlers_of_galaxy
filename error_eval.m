function err = error_eval(X_k,x_vec,s,type)

% Evaluate the error function 
% Input
% X_k: Query state, e.g. R_k,Theta_k 
% x_vec: Vector of the states of settled stars, e.g. r_i or theta_i
% s: Threshold state position or angle
% type: 1 for r, 2 for theta
%
% Output
% err value

    
err =0;

for i=1:length(X_k)
    f = f_eval(X_k(i), x_vec,s);
    g = g_eval(X_k(i), type);
    err = err + (f/g-1)^2;
end
