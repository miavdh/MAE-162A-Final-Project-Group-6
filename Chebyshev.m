% ChebyShev Linkage At Rightmost Deflection
function [min2, min3, min4, h, dfl] = Chebyshev(t2, r1, r2, r3, r4, t1)

% Define variables
theta3 = 0;
theta4 = 0;
eps = 1e-6;

% Number of guesses between 0 and 2pi
numguess = 8;

% Set minimum angle guess for 3 and 4
min2 = 360;
min3 = 360;
min4 = 360;

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
        J = [-r3 * sin(n_1(1,1)), -r4 * sin(n_1(2,1)); r3 * cos(n_1(1,1)), r4 * cos(n_1(2,1))]; % partial of r2 w.r.t y
        f(1,1) = r2 * cos(t2) + r4 * cos(n_1(2,1)) + ...
                 r3 * cos(n_1(1,1)) - r1 * cos(t1);
        f(2,1) = r2 * sin(t2) + r4 * sin(n_1(2,1)) + ...
                 r3 * sin(n_1(1,1));
        
        % Eliminate singular matrices that will not return valid solutions
        if cond(J) > 1e10
            continue
        end
        
        % Newton-Raphson's
        while ((difx > eps || dify > eps) && iter < 1000)            
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
            
            J(1,1) = -r3 * sin(n_1(1,1)); % partial of f1 w.r.t x
            J(1,2) = -r4 * sin(n_1(2,1)); % partial of f2 w.r.t x
            J(2,1) = r3 * cos(n_1(1,1)); % partial of f1 w.r.t y
            J(2,2) = r4 * cos(n_1(2,1)); % partial of r2 w.r.t y

            iter = iter + 1;
        end
            
        % Determine theta2
        % theta2 = asin((4 * A - r3*sin(n(1,1))/2)/r2);

        % Determine value of angles within 0 to 2pi
        true2 = 180 * wrapTo2Pi(t2) / pi();
        true3 = 180 * wrapTo2Pi(n_1(1,1)) / pi();
        true4 = 180 * wrapTo2Pi(n_1(2,1)) / pi();

        % Find lowest value of theta 3 (angle should be acute)
        if (true3 < min3)
            min2 = true2;
            min3 = true3;
            min4 = true4;
        end
    end
end

% Determine maximum hzn and vert deflection 
dfl = r2 * cos(t2) + r3 * cosd(min3) / 2 - r1 / 2;
h = r2 * sin(t2) + r3 * sind(min3) / 2;

end