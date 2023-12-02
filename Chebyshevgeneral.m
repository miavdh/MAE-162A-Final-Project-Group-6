% ChebyShev Linkage At Any Deflection
function [t3, t4, h, dfl] = Chebyshevgeneral(t2, A, t3init, t4init)
% Clear cache
% clear
% close all
% clc

% Test variables
r1 = 4*A;
r2 = 5*A;
r3 = 2*A;
r4 = r2;
% t2 = 36.899*pi()/180;
        
% Initialize variables
difx = 1;
dify = 1;
iter = 0;
min3 = 0;
min4 = 0;
% theta3 = 2 * i * pi() / numguess;
% theta4 = 2 * j * pi() / numguess;
% n_1 = [t3init; t4init];

% Initialize J and f
% J = [-r3 * sin(n_1(1,1)), -r4 * sin(n_1(2,1)); r3 * cos(n_1(1,1)), r4 * cos(n_1(2,1))];
% f(1,1) = r2 * cos(t2) + r4 * cos(t4init) + ...
%          r3 * cos(t3init) - r1;
% f(2,1) = r2 * sin(t2) + r4 * sin(t4init) + ...
%          r3 * sin(t3init);

% Number of guesses between 0 and 2pi
numguess = 8;

% Set range for angle of 3
t3min = pi/2;
t3max = 7*pi/4;

% Run through range of guesses between 0 and 2pi to get complete picture
% of solutions
for i = 1:numguess
    for j = 1:numguess

        % Initialize variables
        difx = 1;
        dify = 1;
        iter = 0;
        theta3 = 2 * i * pi() / numguess;
        theta4 = 2 * j * pi() / numguess;
        n_1 = [theta3; theta4];

        % Initialize J and f
        J = [-r3 * sin(n_1(1,1)), -r4 * sin(n_1(2,1)); r3 * cos(n_1(1,1)), r4 * cos(n_1(2,1))];
        f(1,1) = r2 * cos(t2) + r4 * cos(n_1(2,1)) + ...
                 r3 * cos(n_1(1,1)) - r1;
        f(2,1) = r2 * sin(t2) + r4 * sin(n_1(2,1)) + ...
                 r3 * sin(n_1(1,1));
        
        % Eliminate singular matrices that will not return valid solutions
        if cond(J) > 1e10
            continue
        end
        

        % Newton-Raphson's
        while ((difx > eps || dify > eps) && iter < 1000)
            % Calculate inverse of J
            % Jinv = inv(J);
            
            % Determine value of xn based on xn-1
            n = -J\f + n_1;
            
            % Calculate difference between xn and xn-1
            difx = abs(n(1,1) - n_1(1,1));
            dify = abs(n(2,1) - n_1(2,1));
            
            % New xn-1 is current xn
            n_1 = n;
        
            % Print value of xn
            %fprintf('Iteration: %1.0f \n New xn-1: %1.6f, %1.6f \n', iter, n(1,1), n(2,1));
        
            % Recalculating values of f1, f2, and J
            f(1,1) = r2 * cos(t2) + r4 * cos(n_1(2,1)) + ...
                     r3 * cos(n_1(1,1)) - r1;
            f(2,1) = r2 * sin(t2) + r4 * sin(n_1(2,1)) + ...
                     r3 * sin(n_1(1,1));
            
            J = [-r3 * sin(n_1(1,1)), -r4 * sin(n_1(2,1)); r3 * cos(n_1(1,1)), r4 * cos(n_1(2,1))]; 
        
            iter = iter + 1;
        end

        t3 = wrapTo2Pi(n_1(1,1));
        t4 = wrapTo2Pi(n_1(2,1));

        % Find lowest value of theta 3 (angle should be acute)
        if (t3 < t3max && t3 > t3min)
            min3 = t3;
            min4 = t4;
        end

        %fprintf("%1.4f    %1.4f    %1.5f    %1.5f    %1.5f\n", theta3, theta4, true3, true4, true2);
    end
end

t3 = wrapTo2Pi(min3);
t4 = wrapTo2Pi(min4);

% Determine maximum deflection and 
h = r2 * sin(t2) + r3 * sin(t3) / 2;

end