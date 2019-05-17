function g_val = g_eval(X_k,type)

% Evaluate the f function 
% Input
% X_k: Query state, e.g. R_k,Theta_k 
%
% type: 1 for r, 2 for theta
%
% Output
% g function

global R_min R_max



if type == 1
    r = X_k;
    % Calculate alpha(r)
    if r==2 % kpc
        alpha = 0.5833; 
    elseif r==32 % kpc
        alpha = 0.4948;
    else
        alpha = 1;
    end
    
    g_val = alpha * 2 * r / (R_max^2-R_min^2);
    
elseif type == 2 
    theta = X_k;
    % Calculate beta(theta)
    if theta == -pi
        beta = 0.5;
    elseif theta == pi
        beta = 0.5;
    else 
        beta = 1;
    end
        
    g_val = (1/(2*pi)) * beta;
end
    
    