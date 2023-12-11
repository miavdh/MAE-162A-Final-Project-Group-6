function [t5p, t5pp, t6p, t6pp] = Loop2(r2, r26, r15, r5, r6, t2, t3, t4, t5, t6, t3p, t4p, t3pp, t4pp)

% Calculate Jacobian
J = [-r5*sin(t5), -r6*sin(t6);
     r5*cos(t5), r6*cos(t6)];


% Solve for first order kinematic coefficients
KCv = J\[(r2*sin(t2)+r26*t3p*sin(t3)+r15*t4p*sin(t4)); 
    (-r2*cos(t2)-r26*t3p*cos(t3)-r15*t4p*cos(t4))];
% isolate first order kinematic coefficients
t5p = KCv(1,1);
t6p = KCv(2,1);

% Solve for second order kinematic coefficients
KCa = J\[(r2*cos(t2)+r26*t3pp*sin(t3)+r26*t3p^2*cos(t3)+r6*t6p^2*cos(t6)+...
    r5*t5p^2*cos(t5)+r15*t4pp*sin(t4)+r15*t4p^2*cos(t4));
    (r2*sin(t2)-r26*t3pp*cos(t3)+r26*t3p^2*sin(t3)+r6*t6p^2*sin(t6)+...
    r5*t5p^2*sin(t5)-r15*t4pp*cos(t4)+r15*t4p^2*sin(t4));];
% isolate second order kinematic coefficients
t5pp = KCa(1,1);
t6pp = KCa(2,1);
