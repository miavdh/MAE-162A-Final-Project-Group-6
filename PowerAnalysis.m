function torque = PowerAnalysis(r2, r3, r4, r5, r6, h6, ...
     angle2, angle3, angle4, angle5, angle6, theta31, theta32, ...
     theta41, theta42, theta51, theta52, theta61, theta62, omega2, alpha2)

% Define link lengths - convert from mm to m
r2 = r2/1000;
r3 = r3/1000;
r4 = r4/1000;
r5 = r5/1000;
r6 = r6/1000;
h6 = h6/1000;

% Define inputs
% Calculate mass
rho = 2698.9; % kg/m^3 (density)
% rho = rho / 1000^3; % kg/mm^3
dL = 30/1000; % thickness of link
wL = 30/1000; % width of link

% Actual lengths
r2a = r2 + wL;
r3a = r3 + wL;
r4a = r4 + wL;
r5a = r5 + wL;
r6a = r6 + wL;

% Mass in kg
m2 = rho*r2a*wL*dL;
m3 = rho*r3a*wL*dL;
m4 = rho*r4a*wL*dL;
m5 = rho*r5a*wL*dL;
m6 = rho*r6a*wL*dL + 2*rho*(h6-wL)*wL*dL;

% Calculate moments of inertia
i2 = (1/12)*m2*((wL^2)+(r2)^2);
i3 = (1/12)*m3*((wL^2)+(r3)^2);
i4 = (1/12)*m4*((wL^2)+(r4)^2);
i5 = (1/12)*m5*((wL^2)+(r5)^2);
i6 = (1/12)*m6*((wL^2)+(r6)^2);

%% Repeat for 3 periods of oscillation

% Find KC of center of mass of link 2, 3, 4, 5, 6
    % key: xCoMrj = x position of center of mass of link j
    % yCoMrj = y position of center of mass of link j
    % x1CoMrj = first-order kinematic coeff. of x position of CoM of link j
    % y1CoMrj = first-order kinematic coeff. of y position of CoM of link j
    % x2CoMrj = second-order kinematic coeff. of x position of CoM of link j
    % y2CoMrj = second-order kinematic coeff. of y position of CoM of link j
    % CoMaj = angle of center of mass for link j


xCoMr2 = 0.5*r2*cos(angle2);
yCoMr2 = 0.5*r2*sin(angle2);
x1CoMr2 = -0.5*r2*sin(angle2);
y1CoMr2 = 0.5*r2*cos(angle2);
x2CoMr2 = -xCoMr2;
y2CoMr2 = -yCoMr2;
% CoMa2 = angle2;

xCoMr3 = (r2*cos(angle2))+(0.5*r3*cos(angle3));
yCoMr3 = (r2*sin(angle2))+(0.5*r3*sin(angle3));
x1CoMr3 = (-r2*sin(angle2))-(0.5*theta31*r3*sin(angle3));
y1CoMr3 = (r2*cos(angle2))+(0.5*theta31*r3*cos(angle3));
x2CoMr3 = (-r2*cos(angle2))-(0.5*theta32*r3*sin(angle3))-(0.5*(theta31)^2*r3*cos(angle3));
y2CoMr3 = (-r2*sin(angle2))+(0.5*theta32*r3*cos(angle3))-(0.5*(theta31)^2*r3*sin(angle3));
% CoMa3 = atan(yCoMr3/xCoMr3);

xCoMr4 = (r2*cos(angle2))+(r3*cos(angle3))+(0.5*r4*cos(angle4));
yCoMr4 = (r2*sin(angle2))+(r3*sin(angle3))+(0.5*r4*sin(angle4));
x1CoMr4 = (-r2*sin(angle2))-(r3*theta31*sin(angle3))-(0.5*r4*theta41*sin(angle4));
y1CoMr4 = r2*cos(angle2) + r3*theta31*cos(angle3) + 0.5*r4*theta41*cos(angle4);
x2CoMr4 = -r2*cos(angle2)-r3*theta32*sin(angle3)-r3*theta31^2*cos(angle3)-0.5*r4*theta42*sin(angle4)-0.5*theta41^2*r4*cos(angle4);
y2CoMr4 = -r2*sin(angle2)+r3*theta32*cos(angle3)-r3*theta31^2*sin(angle3)+0.5*r4*theta42*cos(angle4)-0.5*theta41^2*r4*sin(angle4);
% CoMa4 = atan((yCoMr4)/(xCoMr4));

xCoMr5 = xCoMr3+(r6*cos(angle6))+(0.5*r5*cos(angle5));
yCoMr5 = yCoMr3+(r6*sin(angle6))+(0.5*r5*sin(angle5));
x1CoMr5 = x1CoMr3-r6*theta61*sin(angle6)-0.5*r5*theta51*sin(angle5);
y1CoMr5 = y1CoMr3+r6*theta61*cos(angle6)+0.5*r5*theta51*cos(angle5);
x2CoMr5 = x2CoMr3-r6*theta62*sin(angle6)-r6*theta61^2*cos(angle6)-0.5*r5*theta52*sin(angle5)-0.5*r5*theta51^2*cos(angle5);
y2CoMr5 = y2CoMr3+r6*theta62*cos(angle6)-r6*theta61^2*sin(angle6)+0.5*r5*theta52*cos(angle5)-0.5*r5*theta51^2*sin(angle5);
% CoMa5 = atan(yCoMr5/xCoMr5);

% Center of mass of link 6 relative to the tail of r6
xCoMr6_t = 202.144/1000; % found from solidworks
yCoMr6_t = 163.522/1000; % found from solidworks
CoMr6 = sqrt(xCoMr6_t^2 + yCoMr6_t^2); % distance
CoMa6_t = atan(yCoMr6_t/xCoMr6_t); % angle

xCoMr6 = xCoMr3 + CoMr6*cos(angle6+CoMa6_t);
yCoMr6 = yCoMr3 + CoMr6*sin(angle6+CoMa6_t);
x1CoMr6 = x1CoMr3-CoMr6*theta61*sin(angle6+CoMa6_t);
y1CoMr6 = y1CoMr3+CoMr6*theta61*cos(angle6+CoMa6_t);
x2CoMr6 = x2CoMr3-CoMr6*theta62*sin(angle6+CoMa6_t)-CoMr6*theta61^2*cos(angle6+CoMa6_t);
y2CoMr6 = y2CoMr3+CoMr6*theta62*cos(angle6+CoMa6_t)-CoMr6*theta61^2*sin(angle6+CoMa6_t);
% CoMa6 = atan(yCoMr6/xCoMr6);


%  Given constant omega
% w = 1 * pi / 180;
% thetamax = thetamax*pi/180;
% thetamin = thetamin*pi/180;

% Calculate the Torque Needed at any Time 
% P = T2 * w2 = sum(Aj * w2 * alpha2) + sum(Bj * w2^2) + sum(mj * g * y1CoMrj * w2)
% Can cancel w2 on both sides
% T2 = sum(Aj * alpha2) + sum(Bj * w2^2) + sum(mj * g * y1CoMrj)

% Calculate the Change in Kinetic Energy dT/dt 
A2 = Aj(m2,x1CoMr2,y1CoMr2, i2, 1);
A3 = Aj(m3,x1CoMr3,y1CoMr3, i3, theta31);
A4 = Aj(m4,x1CoMr4,y1CoMr4, i4, theta41);
A5 = Aj(m5,x1CoMr5,y1CoMr5, i5, theta51);
A6 = Aj(m6,x1CoMr6,y1CoMr6, i6, theta61);

B2 = Bj(m2, x1CoMr2, x2CoMr2, y1CoMr2, y2CoMr2, i2, 1, 0);
B3 = Bj(m3, x1CoMr3, x2CoMr3, y1CoMr3, y2CoMr3, i3, theta31, theta32);
B4 = Bj(m4, x1CoMr4, x2CoMr4, y1CoMr4, y2CoMr4, i4, theta41, theta42);
B5 = Bj(m5, x1CoMr5, x2CoMr5, y1CoMr5, y2CoMr5, i5, theta51, theta52);
B6 = Bj(m6, x1CoMr6, x2CoMr6, y1CoMr6, y2CoMr6, i6, theta61, theta62);

A = (A2 + A3 + A4 + A5 + A6) * alpha2; % divide both sides by omega2
B = (B2 + B3 + B4 + B5 + B6) * omega2^2;

dTdt = A + B;

% Calculate the Change in Gravitational Potential Energy dU/dt

Ugrav2 = m2 * 9.81 * y1CoMr2;
Ugrav3 = m3 * 9.81 * y1CoMr3;
Ugrav4 = m4 * 9.81 * y1CoMr4;
Ugrav5 = m5 * 9.81 * y1CoMr5;
Ugrav6 = m6 * 9.81 * y1CoMr6;
dUdt= Ugrav2 + Ugrav3 + Ugrav4 + Ugrav5 + Ugrav6;

% Calculate Torque at a Time t

torque = dUdt + dTdt;
