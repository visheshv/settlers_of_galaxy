function [value,terminate,direction] = statesEvents(t,y)

kpc2km = 30856775814671900;
myr = 1e6*31557600;
kms2kpcpmyr = myr/kpc2km;

r = y(1:3);
v = y(4:6);

rcap = r/norm(r);
kcap = [0;0;-1];
vcap = cross(kcap,rcap);

k = [0.00287729 0.0023821 -0.0010625 0.000198502 -1.88428e-05 9.70521e-07 -2.70559e-08 3.7516e-10 -1.94316e-12];
rnorm = norm(r);
vmag = kms2kpcpmyr*([rnorm.^0 rnorm.^1 rnorm.^2 rnorm.^3 rnorm.^4 rnorm.^5 rnorm.^6 rnorm.^7 rnorm.^8]*k').^-1; % Star speeds kpc/Myr
v_star_vec = vcap*vmag;

delv =  (v_star_vec - v)*(1/kms2kpcpmyr);
delv_max = 300;

if ((delv_max < (norm(delv))  && t >= 1) || (rnorm<2) || (rnorm>32)) && rem(t,0.5)<1e-12 % check delV for capture 1Myr after the ms ship burn and @ 0.5 Myr rate
    value = 0 ;
else
    value = 1;
end

terminate = 0;   % Stop the integration
direction =  0;