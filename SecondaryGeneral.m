% Secondary Linkage At Any Deflection
function [t5, t6] = SecondaryGeneral(A, r15, r5, r6, t2, t3, t4, t5init, t6init)
% Clear cache
% clear
% close all
% clc

% Test variables
r1 = 4*A;
r2 = 5*A;
r3 = 2*A;
r26 = r3/2;
        
% Initialize variables
difx = 1;
dify = 1;
iter = 0;
% theta3 = 2 * i * pi() / numguess;
% theta4 = 2 * j * pi() / numguess;
n_1 = [t5init; t6init];

% Initialize J and f
J = [-r5*sin(n_1(1,1)), -r6*sin(n_1(2,1)) ; r5*cos(n_1(1,1)), r6*cos(n_1(2,1))];
f(1,1) = r2*cos(t2) + r26*cos(t3) + r6*cos(n_1(2,1)) + r5*cos(n_1(1,1)) + r15*cos(t4) - r1;
f(2,1) = r2*sin(t2) + r26*sin(t3) + r6*sin(n_1(2,1)) + r5*sin(n_1(1,1)) + r15*sin(t4);

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
    f(1,1) = r2*cos(t2) + r26*cos(t3) + r6*cos(n_1(2,1)) + r5*cos(n_1(1,1)) + r15*cos(t4) - r1;
    f(2,1) = r2*sin(t2) + r26*sin(t3) + r6*sin(n_1(2,1)) + r5*sin(n_1(1,1)) + r15*sin(t4);
    
    J = [-r5*sin(n_1(1,1)), -r6*sin(n_1(2,1)) ; r5*cos(n_1(1,1)), r6*cos(n_1(2,1))];

    iter = iter + 1;
end

% Final values of t5 and t6
t5 = wrapToPi(n_1(1,1));
t6 = wrapToPi(n_1(2,1));

end