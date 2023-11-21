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


    
end