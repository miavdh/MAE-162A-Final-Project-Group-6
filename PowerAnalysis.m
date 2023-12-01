% Full Power Analysis of Linkage

% Define link lengths
A = 200;
r1 = 4*A;
r2 = 5*A;
r3 = 2*A;
r4 = r2;
r15 = 499;
r5 = 501;
r6 = 400;
t1 = 0; % r1 is always horizontal to the right

% Define inputs
% Calculate mass
rho = 2698.9; % kg/m^3 (density)
rho = rho / 1000^3; % kg/mm^3
dL = 30; % thickness of link
wL = 30; % width of link

m1 = rho*(r1+wL)*wL*dL; % kg
m2 = rho*(r2+wL)*wL*dL;
m3 = rho*(r3+wL)*wL*dL;
m4 = rho*(r4+wL)*wL*dL;
m5 = rho*(r5+wL)*wL*dL;

% Calculate moments of inertia

%% Repeat for 3 periods of oscillation

% Position Analysis of Vector Loop 1

% Position Analysis of Vector Loop 2
 
% Velocity Analysis of Vector Loop 1

% Velocity Analysis of Vector Loop 2

% Acceleration Analysis of Vector Loop 1

% Acceleration Analysis of Vector Loop 2

% Find KC of center of mass of link 2, 3, 4, 5, 6