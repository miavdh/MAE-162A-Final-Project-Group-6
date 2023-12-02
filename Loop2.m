function [t5p, t5pp, t6p, t6pp] = Loop2(A, r15, r5, r6, t2, t3, t4, t5, t6)

J = [-r5 * sin(n_2(1,1)), -r6 * sin(n_2(2,1));...
     r5 * cos(n_2(1,1)), r6 * cos(n_2(2,1))]; % partial of r2 w.r.t y


% Second order kinematic coefficients
KC = J/[r2*sin(t2) + KC3_1*r26*sin(t3) + KC4_1*r15*sin(t4) ; 
                -r2*cos(t2) - KC3_1*r26*cos(t3) - KC4_1*r15*cos(t4)];
