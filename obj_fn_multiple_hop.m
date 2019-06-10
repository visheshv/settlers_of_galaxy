function cost = obj_fn_multiple_hop(x0,star_x0,xt,tof)

% Input 
% x0 = design initial state of the settlership as a colmn vector
% star_x0 = initial state of the star as a colmn vector
% xt = final state of the target star as a column vector
% tof = total time of flight of the transfer

% Output 
% cost: total dV (kpc/myr)

tspan1 = linspace(0,x0(13),100);
opts = odeset('AbsTol',1e-5);
[~,states1]=ode45(@state_dynamics,tspan1,x0(1:6,1),opts);

tspan2 = linspace(0,x0(14),100);
[~,states2]=ode45(@state_dynamics,tspan2,x0(7:12,1),opts);

%cost = norm(x0(4:6) - star_x0(4:6)) + norm(states1(end,4:6)'-x0(10:12))  + norm(states2(end,4:6)' - xt(4:6));
cost = x0(13)+x0(14);
dep = norm(x0(4:6) - star_x0(4:6));
inter1 =  norm(states1(end,4:6)'-x0(10:12));
arr = norm(states2(end,4:6)' - xt(4:6));

end