function [t3p, t4p, t3pp, t4pp] = Loop1(r2, r3, r4, t2, t3, t4)
% Velocity analysis of vector loop one

%Calculate the jacobian
J = [-r3*sin(t3), -r4*sin(t4); r3*cos(t3), r4*cos(t4)];

%solve for the first order KC values
KCv = J\[r2*sin(t2) ; -r2*cos(t2)];

% isolate first order kinematic coefficients
t3p = KCv(1,1);
t4p = KCv(2,1);

% solve for acceleration
KCa = J\[r2*cos(t2) + t3p^2*r3*cos(t3) + t4p^2*r4*cos(t4); 
    r2*sin(t2) + t3p^2*r3*sin(t3) + t4p^2*r4*sin(t4)];

% isolate second order kinematic coefficients
t3pp = KCa(1,1);
t4pp = KCa(2,1);
