function [idx, t_departure,t_arrival,is_bad_solution,x0_departure] = find_grid_stars_strategy3(x0,n,star_ID, t_min_departure,query_type,current_star_ID,J_merit,delV_max,delV_used,x,y,z,vx,vy,vz,star_data,star_target_grids)

% Input
% starID: current settled Star IDs,
% t: minimum departure time 
% n: number of target stars
% query_type: 1 for mothership, 2 for settler ship

% current_star_ID: find the ID for the current star to be departed ( !only for settler ships! )
% star_target_grids

% Output
% t_min_departure
% tof
% Target IDs for the 'n' best stars
% is_bad_solution: equals 0 if legitimately good solution found
% x0 state at updated departure conditions

kpc2km = 30856775814671900; 
myr = 1e6*31557600;

% For selecting tof, tof polynomial function versus r (kpc)
poly_n=[1.64965819444119e-08,-1.70149835795047e-06,6.54110993481119e-05,-0.00110272848583409,0.00589360129400973,0.0467622383963570,-0.668482140993974,3.15728858911639,-2.33996016100402];

r_query = norm(x0(1:3));
tof_fit=polyval(poly_n,r_query);    % guess tof in myr
n_fit= (180/16)/tof_fit;            % angular velocity (deg/ myr)

% for
i_tof_margin = round(tof_fit/0.5+1); % number
i_departure = t_min_departure/0.5+1;
i_arrival = min(i_departure+i_tof_margin,181);
t_departure = (i_departure-1) * 0.5;
t_arrival = (i_arrival-1) *0.5;

%% Assign initial positions at t_departure and not t_departure_min
if query_type == 1
    
    r0=x0(1:3);
    v0=x0(4:6);
    x0_departure = [r0;v0];
    
elseif query_type ==2
    
    r0= [x(current_star_ID+1,i_departure),y(current_star_ID+1,i_departure),z(current_star_ID+1,i_departure)]';
    v0= [vx(current_star_ID+1,i_departure),vy(current_star_ID+1,i_departure),vz(current_star_ID+1,i_departure)]';
    x0_departure = [r0;v0];
    
end

%% Find the current grid 

% What could have been the theta at t=0
current_theta= atan3(r0(2),r0(1));
d_theta = n_fit * t_departure * pi/ 180;
init_theta= current_theta + d_theta;
init_theta=mod(init_theta,2*pi);

if init_theta>pi
    init_theta = -pi+(init_theta-pi);   % init_theta represented as an atan2 equivalent expression
end

j = ceil(norm(r0)-2);

if init_theta == -pi
    k =1;
else
    k = ceil( (init_theta +pi) *(16 /pi) );
end

id = star_target_grids(j,k);

% Control k shift to decrease DV transfers
k_shift =0;
tof_extra = 5;

% grid 1
tof_n1= (k_shift+1)* polyval(poly_n,(floor(norm(r0))+1.5)) +tof_extra;    % guess tof in myr
n_grid_1= (180/16)/tof_n1;                       % angular velocity (deg/ myr)
k_1 = mod(k-ceil( abs((n_fit-n_grid_1) * (t_departure + tof_n1) *(pi/180) /(pi/16))),32);
j_1 = j+1;

if k_1 == 0
    k_1 =32;
end

id_1 = star_target_grids(j_1,k_1);


% grid 2
tof_n2= (k_shift+2) * polyval(poly_n,(floor(norm(r0))+1.5))+tof_extra;    % guess tof in myr
k_2 = mod(k_1-1,32);
j_2 = j+1;

if k_2 == 0
    k_2=32;
end

id_2 = star_target_grids(j_2,k_2);


% grid 3
tof_n3=polyval(poly_n,(floor(norm(r0))+0.5))+tof_extra;    % guess tof in myr;    % guess tof in myr
k_3 = mod(k-1,32);
j_3 = j;

if k_3 == 0
    k_3=32;
end

id_3 = star_target_grids(j_3,k_3);

% grid 4
tof_n4= 2 * polyval(poly_n,(floor(norm(r0))+0.5))+tof_extra;    % guess tof in myr;    % guess tof in myr
k_4 = mod(k-2,32);
j_4 = j;

if k_4 == 0
    k_4=32;
end

id_4 = star_target_grids(j_4,k_4);

% grid 5
tof_n5=polyval(poly_n,(floor(norm(r0))-0.5))+tof_extra;    % guess tof in myr
n_grid_5= (180/16)/tof_n5;                       % angular velocity (deg/ myr)
k_5 = mod(k-floor(abs((n_fit-n_grid_5) * t_departure*(pi/180)  /(pi/16))),32);

if k_5 == 0
    k_5=32;
end

j_5 = j+1;
id_5 = star_target_grids(j_5,k_5);

% grid 6
tof_n6=polyval(poly_n,(floor(norm(r0))-0.5))+tof_extra;    % guess tof in myr
n_grid_6= (180/16)/tof_n6;                       % angular velocity (deg/ myr)
k_6 = mod(k-floor(abs((n_fit-n_grid_6) * t_departure*(pi/180)  /(pi/8))),32);

if k_6 == 0
    k_6=32;
end

j_6 = j+1;
id_6 = star_target_grids(j_6,k_6);


for i=1:6
    eval(['tof_n' num2str(i) '=' 'tof_n' num2str(i) '-mod(tof_n' num2str(i) ',0.5)+0.5;']) % Make all tofs multiples of 0.5
end

star_id_settled = star_ID(star_ID~=0); % IDs for settled stars

cond1 = length(star_id_settled) - length(setdiff(star_id_settled,id_1-1));
cond2 = length(star_id_settled) - length(setdiff(star_id_settled,id_2-1));
cond3 = length(star_id_settled) - length(setdiff(star_id_settled,id_3-1));
cond4 = length(star_id_settled) - length(setdiff(star_id_settled,id_4-1));
cond5 = length(star_id_settled) - length(setdiff(star_id_settled,id_5-1));
cond6 = length(star_id_settled) - length(setdiff(star_id_settled,id_6-1));

decision_matrix=[ 0 0 0 0 2 3 1;
                  0 0 0 1 3 2 1;
                  1 0 0 0 3 2 4;
                  0 1 0 0 3 4 1;
                  0 0 1 0 4 2 1;
                  1 1 0 0 3 4 1;
                  0 1 1 0 4 1 3;
                  0 0 1 1 2 1 3;
                  1 0 1 0 4 2 3;
                  1 0 0 1 3 2 4;
                  0 1 0 1 3 1 4;
                  0 1 1 1 1 3 4;
                  1 0 1 1 2 3 4;
                  1 1 0 1 3 4 2;
                  1 1 1 0 4 3 2;
                  1 1 1 1 3 4 2];
              
query_cond = [cond1 cond2 cond3 cond4];
num_zeros = length(query_cond(query_cond==0));
is_bad_solution= max(3-num_zeros,0);
idx_row= find(((decision_matrix(:, 1) == query_cond(1) & decision_matrix(:, 2) == query_cond(2) & decision_matrix(:,3) == query_cond(3)) & decision_matrix(:,4) == query_cond(4))==1);
eval(['idx =[id_' num2str(decision_matrix(idx_row,5)) ',id_' num2str(decision_matrix(idx_row,6)) ',id_' num2str(decision_matrix(idx_row,7)) '];'])
eval(['t_arrival =[t_departure+tof_n' num2str(decision_matrix(idx_row,5)) ',t_departure+tof_n' num2str(decision_matrix(idx_row,6)) ',t_departure+tof_n' num2str(decision_matrix(idx_row,7)) '];'])
eval(['is_bad_solution=' num2str(is_bad_solution) ';']);         
idx =idx-1;

% %% 
% if cond1 == 0 && cond2 == 0 && cond3 == 0
%    idx = [id_3 id_2 id_1] - 1;
%    t_arrival = [ t_departure+tof_n3 t_departure+tof_n2 t_departure+tof_n1];
%    is_bad_solution=0;
% elseif cond1 == 0 && cond2 == 0 && cond3 ~= 0
%    idx = [id_3 id_2 id_1] - 1;
%    t_arrival = [ t_departure+tof_n3 t_departure+tof_n2 t_departure+tof_n1];
%    is_bad_solution=0;
% elseif cond1 ~= 0 && cond2 == 0 && cond3 == 0 && cond4 == 0
%     idx = [id_3 id_4 id_2] - 1;
%     t_arrival = [t_departure+tof_n3 t_departure+tof_n4 t_departure+tof_n2  ];
%     is_bad_solution=0;
% elseif cond1 ~= 0 && cond2 ~= 0 && cond3 == 0 && cond4 == 0
%     idx = [id_3 id_4 id_2] - 1;
%     t_arrival = [t_departure+tof_n3 t_departure+tof_n4 t_departure+tof_n2];
%     is_bad_solution=1;
% elseif cond1 ~= 0 && cond2 ~= 0 && cond3 ~= 0 && cond4 == 0
%     idx = [id_4 id_3 id_2] - 1;
%     t_arrival = [t_departure+tof_n4 t_departure+tof_n3 t_departure+tof_n2];
%     is_bad_solution=2;
% elseif cond1 ~= 0 && cond2 ~= 0 && cond3 ~= 0 && cond4 ~= 0
%     idx = [id_4 id_3 id_2] - 1;
%     t_arrival = [t_departure+tof_n4 t_departure+tof_n3 t_departure+tof_n2];
%     is_bad_solution=3;
% end
%% 
t_arrival(t_arrival>90)=90;

if query_type ==1
    idx=idx(1);
    t_arrival=t_arrival(1);
else
    t_departure = repmat(t_departure,3,1);    
end

end