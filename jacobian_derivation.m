syms x y z vx vy vz kpc2km myr
syms k0 k1 k2 k3 k4 k5 k6 k7 k8
f1= [vx;vy;vz];
r= sqrt(x^2+y^2+z^2);
f2= -(1/r^2) * ((1/ (k0 + k1 * r + k2 * r^2 + k3*r^3 + k4 * r^4 +k5*r^5 +k6*r^6 + k7*r^7 +k8* r^8))*myr/kpc2km)^2 * [x;y;z] ;
F= [f1;f2];
X= [x;y;z;vx;vy;vz];
J=jacobian(F,X)

%% derivation of velocity gradient 

syms r kpc2km myr
syms k0 k1 k2 k3 k4 k5 k6 k7 k8

v= (1/ (k0 + k1 * r + k2 * r^2 + k3*r^3 + k4 * r^4 +k5*r^5 +k6*r^6 + k7*r^7 +k8* r^8))*myr/kpc2km;
J= jacobian(v,r)