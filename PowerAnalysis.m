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

% Define Kinematic Coefficients
theta31 = 3;
theta32 = 3;
theta41 = 3:
theta42 = 3;
theta51 = 3;
theta52 = 3;
theta61 = 3;
theta62 = 3;

% Define Angles(radians)
angle1 = 0;
angle2 = 30;
angle3 = 30;
angle4 = 30;
angle5 = 30;
angle6 = 30;


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

xCoMr2 = 0.5*r2*cos(angle2);
yCoMr2 = 0.5*r2*sin(angle2);
x1CoMr2 = -0.5*r2*sin(angle2);
y1CoMr2 = 0.5*r2*cos(angle2);
x2CoMr2 = -xCoMr2;
y2CoMr2 = -yComr2;
CoMa2 = angle2;

xCoMr3 = (r2*cos(angle2))+(0.5*r3*cos(angle3));
yCoMr3 = (r2*sin(angle2))+(0.5*r3*sin(angle3));
x1CoMr3 = (-r2*sin(angle2))-(0.5*theta31*r3*sin(angle3));
y1CoMr3 = (r2*cos(angle2))-(0.5*theta31*r3*cos(angle3));
x2CoMr3 = (-r2*cos(angle2))-(0.5*theta32*r3*sin(angle3))-(0.5*(theta31)^2*r3*cos(angle3));
y2CoMr3 = (-r2*sin(angle2))+(0.5*theta32*r3*cos(angle3))-(0.5*(theta31)^2*r3*sin(angle3));
CoMa3 = atan(yComr3/xCoMr3);


xCoMr4 = (r1*cos(angle1))-(0.5*r4*cos(angle4));
yCoMr4 = (r1*sin(angle1))+(0.5*r4*sin(angle4));
x1CoMr4 = (0.5*r4*theta41*sin(angle4));
y1CoMr4 = (-0.5*r4*theta41*cos(angle4));
x2CoMr4 = (0.5*theta42*r4*sin(angle4))+(0.5*(theta41)^2*r4*cos(angle4));
y2CoMr4 = (-0.5*theta42*r4*cos(angle4))+(0.5*(theta41)^2*r4*sin(angle4));
CoMa4 = atan((yComr4)/(xCoMr4));

xCoMr5 = (r1*cos(angle1))-(r15*cos(angle4))-(0.5*r5*cos(angle5));
yCoMr5 = (r1*sin(angle1))-(r15*sin(angle4))-(0.5*r5*sin(angle5));
x1CoMr5 = (r15*theta41*sin(angle4))+(0.5*r5*theta51*sin(angle5));
y1CoMr5 = (-r15*theta41*cos(angle4))-(0.5*r5*theta51*cos(angle5));
x2CoMr5 = (r15*theta42*sin(angle4))+(r15*(theta41)^2*cos(angle4))+(0.5*r5*theta52*sin(angle5))+(0.5*r5*(theta51)^2*cos(angle5));
y2CoMr5 = (-r15*theta42*cos(angle4))+(r15*(theta41)^2*sin(angle4))-(0.5*r5*theta52*cos(angle5))+(0.5*r5*(theta51)^2*sin(angle5));
CoMa5 = atan(yComr5/xCoMr5);


xCoMr6 = (r1*cos(angle1))-(r15*cos(angle4))-(r5*cos(angle5))-(0.5*r6*cos(angle6))+(r66cos(angle6+pi));
yCoMr6 = (r1*sin(angle1))-(r15*sin(angle4))-(r5*sin(angle5))-(0.5*r6*sin(angle6))+(r66sin(angle6+pi));
x1CoMr5 = (r15*theta41*sin(angle4))+(r5*theta51*sin(angle5))+(0.5*r6*theta61*sin(angle6))-(r66*theta61*sin(angle6+pi));
y1CoMr6 = (-r15*theta41*cos(angle4))-(r5*theta51*cos(angle5))-(0.5*r6*theta61*cos(angle6))+(r66*theta61*cos(angle6+pi));
x2CoMr6 = (r15*theta42*sin(angle4))+(r15*(theta41)^2*cos(angle4))+(r5*theta52*sin(angle5))+(r5*(theta51)^2*cos(angle5))+(0.5*r6*theta62*sin(angle6))+(0.5*r6*(theta61)^2*cos(angle6))-(r66*theta62*sin(angle6+pi))-(r66*(theta61)^2*cos(angle6+pi));
y2CoMr6 = (-r15*theta42*cos(angle4))+(r15*(theta41)^2*sin(angle4))-(r5*theta52*cos(angle5))+(r5*(theta51)^2*sin(angle5))-(0.5*r6*theta62*cos(angle6))+(0.5*r6*(theta61)^2*sin(angle6))+(r66*theta62*cos(angle6+pi))-(r66*(theta61)^2*sin(angle6+pi));
CoMa6 = atan(yComr6/xCoMr6);

