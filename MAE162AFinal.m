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

clear
close all
clc

% Set fixed values
A = 202; %mm
r1 = 4*A; %mm
r2 = 5*A;
r3 = 2*A;
r4 = r2;
t1 = 0; %rad
r26 = r3/2;

% First Vector Loop
% Run ChebyShev for right-most deflection angle
[t2deg, t3deg, t4deg, del3, maxTrav] = Chebyshev_max(r1, r2, r3, r4, t1, A);

% Determine maximum travel
maxTrav = abs(maxTrav) + 2*A;

% Convert back to radians for further calculations
t2 = t2deg*pi()/180;
t3 = t3deg*pi()/180;
t4 = t4deg*pi()/180;

%% Second Vector Loop
% Set leftmost deflection angles
t2l = t2;
t3l = t3;
t4l = t4;
t15l = t4l;
t26l = t3l;

% Set rightmost deflection angles
t2r = acos(4/5);
t3r = pi/2;
t4r = pi/2;
t15r = t4r;
t26r = t3r;

% Iterating r5, r6, r15
%r6 = 408; % decided based on geometric constraints (fit in rectangle)
r15min = 0;
r15max = r2;
deltar15 = 1;

% Set initial values for optimal position (least movement in center of
% mass)
CoMmin = 10;
r5right = zeros(1,1);
r15right = zeros(1,1);
t5right = zeros(1,1);
t6right = zeros(1,1);
r6right = zeros(1,1);
CoMright = zeros(1,1);

%test
% r5 = 660;
% r15 = 492;
% [t5, t6] = NR(r1, r2, r3, r4, r5, r6, r26, r15, ...
%     t1, t2r, t3r, t4r, t15r, t26r)
% t5*180/pi()
% t6*180/pi()

%% Check rightmost 
LR = false; %true is left, false is right
t6max = pi/4;
for r15 = r15min:deltar15:r15max
    % Calculate dimensions of line going from center of r3 to endpt. of r15
    % rightLengthx = abs(r1 - r15*cos(t15r) - r2*cos(t2r) - r3 * cos(t3r)/2);% length of triangle to determine range of r5
    % rightLengthy = abs(- r2*sin(t2r) - r3*sin(t3r)/2 - r15*sin(t15r));
    % rightLength = sqrt(rightLengthx^2+rightLengthy^2);
    % rightAngle = 2*pi - (3*pi/2 + atan(rightLengthx/rightLengthy));

    leftLengthx = abs(r1 - r15*cos(t15l) - r2*cos(t2l) - r3 * cos(t3l)/2);
    leftLengthy = abs(- r2*sin(t2l) - r3*sin(t3l)/2 - r15*sin(t15l));
    leftLength = sqrt(leftLengthx^2 + leftLengthy^2);
    leftAngle = atan(leftLengthy/leftLengthx);

    % r6 = (rightLength^2 - leftLength^2)/...
    %      (2*(rightLength*cos(rightAngle) - leftLength*cos(leftAngle)));
    r6 = floor(((4*A-r15)^2 - leftLength^2)/...
         (2*((4*A-r15)*cos(pi/2) - leftLength*cos(leftAngle))));
      
    % Determine range of r5
    phimin = pi/2 - 2*pi/180;
    phimax = pi/2 + 2*pi/180;
    % r5min = sqrt(rightLength^2 + r6^2 - 2*rightLength*r6*cos(phimin));
    % r5max = sqrt(rightLength^2 + r6^2 - 2*rightLength*r6*cos(phimax));
    r5min = sqrt((4*A-r15)^2 + r6^2 - 2*(4*A-r15)*r6*cos(phimin));
    r5max = sqrt((4*A-r15)^2 + r6^2 - 2*(4*A-r15)*r6*cos(phimax));
    r5min = floor(r5min); % make sure lengths are integers
    r5max = ceil(r5max);
    deltar5 = 1;
    for r5 = r5min:deltar5:r5max
        [t5, t6, CoM] = NR(r1, r2, r3, r4, r5, r6, r26, r15, ...
        t1, t2r, t3r, t4r, t15r, LR);
        t6;
        % CoM test   
        if (abs(CoM) < 2 && abs(t6) < t6max)
            if (r5right(1,1) == 0) 
                r5right(1,1) = r5;
                r15right(1,1) = r15;
                t5right(1,1) = t5;
                t6right(1,1) = t6;
                r6right(1,1) = r6;
                CoMright(1,1) = CoM;
            else
                r5right = [r5right; r5];
                r15right = [r15right; r15];
                r6right = [r6right; r6];
                CoMright = [CoMright; CoM];
                % t5right = [t5right; t5];
                % t6right = [t6right; t6];
            end
        end
    end
end

%% Leftmost deflection test
sizeRight = size(r15right);
iterLeft = sizeRight(1,1);

% initialize
r15left = zeros(1,1);
r5left = zeros(1,1);
r6left = zeros(1,1);
t6left = zeros(1,1);

for i = 1:iterLeft
    LR = true;
    r15 = r15right(i,1);
    r5 = r5right(i,1);
    r6 = r6right(i,1);
    [t5, t6, CoM] = NR(r1, r2, r3, r4, r5, r6, r26, r15, ...
                  t1, t2l, t3l, t4l, t15l, LR);
        if (abs(CoM) < 2 && abs(t6) < pi/180)
            if (r5left(1,1) == 0) 
                r5left(1,1) = r5;
                r15left(1,1) = r15;
                r6left(1,1) = r6;
                % t5left(1,1) = t5;
                t6left(1,1) = t6;
            else
                r5left = [r5left; r5];
                r15left = [r15left; r15];
                r6left = [r6left; r6];
                % t5left = [t5left; t5];
                t6left = [t6left; t6];
            end
        end
end

%% Centermost deflection test
sizeLeft = size(r15left);
iterCenter = sizeLeft(1,1);

% initialize
r15center = zeros(1,1);
r5center = zeros(1,1);
r6center = zeros(1,1);
t6center = zeros(1,1);

% Center angles
t2c = asin(4/5);
t3c = pi;
t4c = 2*pi - asin(4/5);
t15c = t4c;
t26c = t3c;

for i = 1:iterCenter
    LR = true; 
    r15 = r15left(i,1);
    r5 = r5left(i,1);
    r6 = r6left(i,1);
    [t5, t6, CoM] = NR(r1, r2, r3, r4, r5, r6, r26, r15, ...
                  t1, t2c, t3c, t4c, t15c, LR);
        if (abs(CoM) < 1 && abs(t6) < pi/180)
            if (r5center(1,1) == 0) 
                r5center(1,1) = r5;
                r15center(1,1) = r15;
                r6center(1,1) = r6;
                % t5left(1,1) = t5;
                t6center(1,1) = t6;
            else
                r5center = [r5center; r5];
                r15center = [r15center; r15];
                r6center = [r6center, r6];
                % t5left = [t5left; t5];
                t6center = [t6center; t6];
            end
        end
end

%% Test different deflection values to see if table moves up and down


% CoMmin
% r5best
% r15best
% t5best
% t6best

% Newton-Raphson's to solve for angles
% function [t5, t6] = NR(r1, r2, r3, r4, r5, r6, r16, r15, ...
%     t1, t2, t3, t4, t15, t16)
% 
% %make D very large beginning so we can start the loop
% d5 = 100;
% d6 = 100;
% D = [d5;d6];
% 
% t5_0 = pi()/180;
% t6_0 = pi()/180;
% 
% %variables for iterations
% n = 0;
% t5_1 = t5_0;
% t6_1 = t6_0;
% 
% while norm(D, 1)>10^(-6);
%     %iteration check
%     n = n + 1; 
% 
%     %solve for the fuctions and jacobian
%     f1 = r2*cos(t2) + r16*cos(t3) + r6*cos(t6_1) + r5*cos(t5_1) + r15*cos(t15) - r1;
%     f2 = r2*sin(t2) + r16*sin(t3) + r6*sin(t6_1) + r5*sin(t5_1) + r15*sin(t4);
%     F = [f1; f2];
% 
%     J = [-r5*sin(t5_1), -r6*sin(t6_1) ; r5*cos(t5_1), r6*cos(t6_1)];
% 
%     % Eliminate singular matrices that will not return valid solutions
%         if cond(J) > 1e10
%             continue
%         end
% 
%     %return values from solved iteration
%     D = -J\F; 
%     d5 = D(1);
%     d6 = D(2);
% 
%     %update the initial parameters
%     t5_1 = t5_1 + d5;
%     t6_1 = t6_1 + d6; 
% end
% 
% t5 = t5_1;
% t6 = t6_1;
% 
% end
