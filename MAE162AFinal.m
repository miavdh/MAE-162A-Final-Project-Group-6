%% MAE 162A Final Project Code
% Names: Ahan Agarwal, Jackson Bullard, Dario Cardenas, Alexi Gill, 
%        Pearl Klassen, Martin Nay, Grace Pelligrino, Rodolfo Ruiz, 
%        Lance Taylor, Mikaela Van de Heetkamp
% Class: MECH&AE 162A
% Instructor: Christopher Matthes
% Due: December 15, 2023 at 10 A.M.
% Goal: Design a mechanism that fits in a 850mm x 1050mm footprint in one
%       configuration. The mechanism must convert an angular input into a
%       (nearly) pure translational output. The translation must be at
%       least 850mm along the x-axis, with less than 4mm of vertical
%       translation of the table's center of mass along the y-axis. In
%       addition, the table cannot tilt more than 1deg above or below the
%       horizontal.

% Set fixed values
A = 212.5; %mm
r1 = 4*A; %mm
r2 = 5*A;
r3 = 2*A;
r4 = r2;
t1 = 0; %rad

% Set leftmost deflection angles
t2l = pi()/2;
t3l = 
t4l = atan(3/4);


% Newton-Raphson's to solve for angles
function [t5, t6] = NR(r1, r2, r3, r4, r5, r6, r16, r15, ...
    t1, t2, t3, t4, t15, t16)

%make D very large beginning so we can start the loop
d5 = 100;
d6 = 100;
D = [d5;d6];

t5_0 = pi()/180;
t6_0 = pi()/180;

%variables for iterations
n = 0;
t5_1 = t5_0;
t6_1 = t6_0;

while norm(D, 1)>10^(-6);
%iteration check
n = n + 1; 

%solve for the fuctions and jacobian
f1 = r2*cos(t2) + r26*cos(t3) + r6*cos(t6) + r5*cos(t5) + r15*cos(t15) - r1;
f2 = r2*sin(t2) + r26*sin(t3) + r6*sin(t6 + r5*sin(t5) + r15*sin(t4);
F = [f1; f2];

J = [-r5*sin(t5), -r6*sin(t6) ; r5*cos(t5), r6*cos(t6)];

%return values from solved iteration
D = -J\F; 
d5 = D(1);
d6 = D(2);

%update the initial parameters
t5_1 = t5_1 + d5;
t6_1 = t6_1 + d6; 
end

t5 = t5_1;
t6 = t6_1;

    
end
