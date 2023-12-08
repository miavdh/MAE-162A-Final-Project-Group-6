% Full Analysis of Linkage
clear
close all
clc

% Define link lengths
A = 202;
r1 = 4*A;
r2 = 5*A;
r3 = 2*A;
r4 = r2;
r15 = 499;
r5 = 507;
r6 = 404;
h6 = 240;
t1 = 0; % r1 is always horizontal to the right
dL = 30; % thickness of link
wL = 30; % width of link
t2min = atan(3/4); 
t2max = 1.662; 

% Set initial values for angles (to begin iterations)
t3init = pi;
t4init = 2*pi-t2min;
t5init = 3*pi/2;
t6init = 0;

% Actual lengths
r1a = r1 + wL;
r2a = r2 + wL;
r3a = r3 + wL;
r4a = r4 + wL;
r5a = r5 + wL;
r6a = r6 + wL;

% Define inputs
% Calculate mass
rho = 2698.9; % kg/m^3 (density)
rho = rho / 1000^3; % kg/mm^3

% Mass in kg
m1 = rho*r1a*wL*dL;
m2 = rho*r2a*wL*dL;
m3 = rho*r3a*wL*dL;
m4 = rho*r4a*wL*dL;
m5 = rho*r5a*wL*dL;
m6 = rho*r6a*wL*dL + 2*rho*(h6-wL)*wL*dL;

% Calculate moments of inertia, kg * mm^4
Ig1 = m1*(wL^2+r1a^2)/12; 
Ig2 = m2*(wL^2+r2a^2)/12;
Ig3 = m3*(wL^2+r3a^2)/12;
Ig4 = m4*(wL^2+r4a^2)/12;
Ig5 = m5*(wL^2+r5a^2)/12;
Ig6 = m6*(wL^2+r6a^2)/12; % PLACEHOLDER

%% Determine full range of positions w.r.t to time, t

% Time limits for position
w = pi/180;
tmin = 0;
tmax = 2*pi/w; % one full period
numStep = 1000;
step = (tmax - tmin)/numStep;

% Set input and output matrices to store values
tFull = zeros(numStep, 1);
t2Full = zeros(numStep, 1); % deflection of 2
t3Full = zeros(numStep, 1);
t4Full = zeros(numStep, 1);
w2Full = zeros(numStep, 1); % ang speed of 2
a2Full = zeros(numStep, 1); % ang accel of 2
dflFull = zeros(numStep, 1); % horizontal deflection of table
CoMFull = zeros(numStep, 1); % vertical deflection of table
t6Full = zeros(numStep, 1);

%indexing variable
i = 0;
for t = tmin:step:tmax
    i = i + 1;
    t2 = ((t2max - t2min)/2) * sin(w*t) + ((t2max + t2min)/2); % angle of 2
    w2 = ((t2max - t2min)/2) * w * cos(w*t); % angular speed of 2
    a2 = -((t2max - t2min)/2) * w^2 * sin(w*t); % angular accel of 2

    % Position Analysis of Vector Loop 1
    % run function
    [t3, t4, h] = Chebyshevgeneral(t2, A, t3init, t4init);

    t3init = t3;
    t4init = t4; % output

    % Position Analysis of Vector Loop 2
    [t5, t6] = SecondaryGeneral(A, r15, r5, r6, t2, t3, t4, t5init, t6init);
    CoM = r2*sin(t2) + r3*sin(t3)/2 + r6*sin(t6) / 2;

    % Set new input values
    t5init = t5;
    t6init = t6;

    % Calculate hzn deflection of COM of link 6
    dfl = r2 * cos(t2) + r3 * cos(t3) / 2 + r6 * cos(t6) / 2 - r1 / 2;

    % Store values
    tFull(i,1) = t;
    t2Full(i,1) = t2; % deflection of 2
    t3Full(i,1) = t3;
    t4Full(i,1) = t4;
    w2Full(i,1) = w2; % ang speed of 2
    a2Full(i,1) = a2; % ang accel of 2
    dflFull(i,1) = dfl; % horizontal deflection of table
    CoMFull(i,1) = CoM; % vertical deflection of table
    t6Full(i,1) = t6;

    %Troubleshooting
    fprintf('t2 = %1.5f, t3 = %1.5f, t4 = %1.5f\n t5 = %1.5f, t6 = %1.5f\n',...
             180*t2/pi,  180*t3/pi,  t4*180/pi,  180*t5/pi,  180*t6/pi);
    fprintf('CoM = %1.5f\n', CoM)
end
vert = max(CoMFull) - min(CoMFull)
t6min = min(t6Full)*180/pi
t6max = max(t6Full)*180/pi
hzn = max(dflFull) - min(dflFull)

figure(1)
plot(t2Full, dflFull)
figure(2)
plot(t2Full, CoMFull)
figure(3)
plot(t2Full, t6Full*180/pi)

%% Repeat for 3 periods of oscillation

%% Determine full range of positions w.r.t to time, t

% Time limits for power
w = pi/180;
tmin = 0;
tmax = 6*pi/w; % three full periods
numStep = 3000;
step = (tmax - tmin)/numStep;

% Set input and output matrices to store values
tPow = zeros(numStep, 1);
t2Pow = zeros(numStep, 1); % deflection of 2
w2Pow = zeros(numStep, 1); % ang speed of 2
a2Pow = zeros(numStep, 1); % ang accel of 2
t3Pow = zeros(numStep, 1);
t4Pow = zeros(numStep, 1);
t5Pow = zeros(numStep, 1);
t6Pow = zeros(numStep, 1);
dflPow = zeros(numStep, 1); % horizontal deflection of table
CoMPow = zeros(numStep, 1); % vertical deflection of table

for t = tmin:step:tmax
    i = i + 1;
    t2 = ((t2max - t2min)/2) * sin(w*t) + ((t2max + t2min)/2); % angle of 2
    w2 = ((t2max - t2min)/2) * w * cos(w*t); % angular speed of 2
    a2 = -((t2max - t2min)/2) * w^2 * sin(w*t); % angular accel of 2

    % Position Analysis of Vector Loop 1
    % run function
    [t3, t4, h] = Chebyshevgeneral(t2, A, t3init, t4init);

    t3init = t3;
    t4init = t4; % output

    % Position Analysis of Vector Loop 2
    [t5, t6] = SecondaryGeneral(A, r15, r5, r6, t2, t3, t4, t5init, t6init);


    % Set new input values
    t5init = t5;
    t6init = t6;

    [torque] = PowerAnalysis(A, tmax, tmin, 
    % Store values
    torque(i,1) = torque;

    
end


% Velocity/Accel Analysis of Vector Loop 1
[t3p, t4p, t3pp, t4pp] = Loop1(A, t2, t3, t4, t, w);

% Velocity/Acceleration Analysis of Vector Loop 2


% Find KC of center of mass of link 2, 3, 4, 5, 6
