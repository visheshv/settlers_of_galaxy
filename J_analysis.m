% Load data
clear all
close all

load star_data

global R_min R_max

R_min = 2; % kpc
R_max = 32; % kpc

R_vec = data(:,2);       % kpc (kiloparsecs)
theta_f_vec = data(:,6); % final polar angle (deg)

figure(1)
histogram(R_vec);
xlabel('position (kpc)');ylabel('num of stars')

% Check the histogram of the positions of stars
figure(2)
histogram(theta_f_vec);
xlabel('final polar angle (deg)');ylabel('num of stars')

% Sort data by radial position: Closest stars first
data_sorted = sortrows(data,2); % Sort ascending by radial position / theta position 2/6 

err_store = zeros(100,1);
error_J_term = zeros(100,1);
count = 0;


for i = 1:1e3:100e3 % Move up by 1000 stars to be settled

count=count+1;
n_stars=25000;
    if i>=(100e3-n_stars)
        data = data_sorted(i:end,:); % filter the first 'i' stars
    else
        data = data_sorted(i:(i+n_stars),:); % filter the first 'i' stars
    end

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

err_store(count) = err_r+err_theta;
error_J_term(count) = i / (1 + 1e-4 * i * err_store(count));

end

figure(3)
plot(err_store)
xlabel('Number of closest stars settled  x 1000')
ylabel('E_r + E_{theta}')

figure(4)
plot(error_J_term)
xlabel('Number of closest stars settled  x 1000')
ylabel(' N / (1 + N x 1e-4 x (E_r + E_{theta}))')


% Analysis on submission delay

tspan = 23; % days;
tdue  = 23:-1:0;
B = 1+ (1/tspan)^4 * tdue .^4;
figure(5)
plot(0:23,B);
xlabel('time submission from tstart');
ylabel('B')

%% check by occupying one star in all grids
close all

j_store=[];
j_merit=0;

for radius = 2:1:31 % Move up by 1000 stars to be settled

    for theta = -pi: pi/16:(pi-pi/16)
        
        cond1=data(:,2)>=radius;
        cond2=data(:,2)<(radius+1);
        cond3=data(:,6)>=theta;
        cond4=data(:,6)<(theta+pi/16);
        
        id = find(cond1 & cond2 & cond3 & cond4);
                
        if length(id)==0
            
            j_merit = j_merit + 0 ;
            j_store = [j_store j_merit];
            
        else
            
            % Lets look at the R_, theta_f data for one arbitrary star
            
            
            [~,id_center] = min((data(id,2)- (radius+0.5)).^2 + (data(id,6)- (theta+pi/32)).^2);
            id = id(id_center);
            
            R_vec = data(id,2);       % kpc (kiloparsecs)
            theta_f_vec = data(id,6); % final polar angle (deg)
            
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
            
            err_store= err_r+err_theta;
            error_J_term= i / (1 + 1e-4 * i * err_store);
            
            j_merit = j_merit + error_J_term ;
            j_store = [j_store j_merit];
        
        end
        
    end
    
end

disp(['Final value of merit function: ' num2str(j_merit)])

figure(6)
plot(j_store)
xlabel('Number of grids covered')
ylabel('J value')





