function [c,ceq] = nonlcon_multiple_hop(x0,star_x0,xt,tof)

tspan1 = linspace(0,x0(13),100);
opts = odeset('AbsTol',1e-5);
[~,states1]=ode45(@state_dynamics,tspan1,x0(1:6,1),opts);

tspan2 = linspace(0,x0(14),100);
[~,states2]=ode45(@state_dynamics,tspan2,x0(7:12,1),opts);

kpc2km = 30856775814671900;myr = 1e6*31557600;kms2kpcpmyr = myr/kpc2km;

c = [norm(x0(4:6) - star_x0(4:6)) - 174.9*kms2kpcpmyr;norm(states1(end,4:6)'-x0(10:12))-174.9*kms2kpcpmyr;norm(states2(end,4:6)' - xt(4:6))-174.9*kms2kpcpmyr];

ceq = [(x0(1:3)-star_x0(1:3)); (states1(end,1:3)'-x0(7:9)); (states2(end,1:3)'-xt(1:3));norm(x0(4:6) - star_x0(4:6)) + norm(states1(end,4:6)'-x0(10:12))  + norm(states2(end,4:6)' - xt(4:6))-399.9*kms2kpcpmyr];

end