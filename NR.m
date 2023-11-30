% Newton-Raphson's to solve for angles
function [t5true, t6true] = NR(r1, r2, r3, r4, r5, r6, r26, r15, t1, t2, t3, t4, t15, t16, LR)
% clc
% clear all
% close

% A = 207; %mm
% r1 = 4*A; %mm
% r2 = 5*A;
% r3 = 2*A;
% r4 = r2;
% t1 = 0; %rad
% 
% % Set rightmost deflection angles
% t2 = atan(3/4);
% t3 = pi()/2;
% t4 = 3*pi()/2;
% r26 = r3/2;
% t15 = t4;
% t26 = t3;
% r5 = 850;

% Iterating r5, r6, r15
% r6 = 371;
% r15 = r4 - sqrt(r5^2 - r6^2);

% Limits of angle 5 based on left or right limit
if (LR == false)
    t5min = pi();
    t5max = 3*pi()/2;
else
    t5min = pi();
    t5max = 2*pi();
end

%% Newton raphsons

% Number of guesses between 0 and 2pi
numguess = 8;

% Set angle guess for 5 and 6
t5true = 2*pi();
t6true = 2*pi();

% Run through range of guesses between 0 and 2pi to get complete picture
% of solutions
for i = 1:numguess
    for j = 1:numguess

        % Initialize variables
        difx = 1;
        dify = 1;
        iter = 0;
        t5 = 2 * i * pi() / numguess;
        t6 = 2 * j * pi() / numguess;
        n_1 = [t5; t6];
        eps = 10^-6;
        
        % Initialize J and f
        J = [-r5*sin(n_1(1,1)), -r6*sin(n_1(2,1)) ; r5*cos(n_1(1,1)), r6*cos(n_1(2,1))];
        f(1,1) = r2*cos(t2) + r26*cos(t3) + r6*cos(n_1(2,1)) + r5*cos(n_1(1,1)) + r15*cos(t15) - r1;
        f(2,1) = r2*sin(t2) + r26*sin(t3) + r6*sin(n_1(2,1)) + r5*sin(n_1(1,1)) + r15*sin(t15);
        
        % Eliminate singular matrices that will not return valid solutions
        if cond(J) > 1e10
            continue
        end
        
        % Newton-Raphson's
        while ((difx > eps || dify > eps) && iter < 1000)
            % Calculate inverse of J
            Jinv = inv(J);
            
            % Determine value of xn based on xn-1
            n = -Jinv * f + n_1;
            
            % Calculate difference between xn and xn-1
            difx = abs(n(1,1) - n_1(1,1));
            dify = abs(n(2,1) - n_1(2,1));
            
            % New xn-1 is current xn
            n_1 = n;
        
            % Print value of xn
            %fprintf('Iteration: %1.0f \n New xn-1: %1.6f, %1.6f \n', iter, n(1,1), n(2,1));
        
            % Recalculating values of f1, f2, and J
            f(1,1) = r2*cos(t2) + r26*cos(t3) + r6*cos(n_1(2,1)) + r5*cos(n_1(1,1)) + r15*cos(t15) - r1;
            f(2,1) = r2*sin(t2) + r26*sin(t3) + r6*sin(n_1(2,1)) + r5*sin(n_1(1,1)) + r15*sin(t15);
            
            J = [-r5*sin(n_1(1,1)), -r6*sin(n_1(2,1)) ; r5*cos(n_1(1,1)), r6*cos(n_1(2,1))];
        
            iter = iter + 1;
        end
        
        % Final values of t5 and t6
        t5 = wrapTo2Pi(n_1(1,1));
        t6 = wrapTo2Pi(n_1(2,1));

        % Test if angles are within expected range and minimize angle of 6
        if (t5 >= t5min && t5 <= t5max && (t6 <= pi/2 || t6 >= 3*pi/2))
            t5true = t5;
            t6true = t6;
        end
    end
end

% function end
end