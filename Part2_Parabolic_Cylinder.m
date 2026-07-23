%% Differential Equations Final Project -- PART 2: Parabolic Cylinder Equation
% shima samadi , mahyas javaheri
% Term: 1404_05 Second Semester | Instructor: Dr. Shameli

clc;
clear;
close all;
fprintf('\n============================================================\n');
fprintf('PART 2: Parabolic Cylinder Equation\n');
fprintf('============================================================\n');
n_vals = [1, 5, 15];
xspan2 = [0, 10];
% 2-2-1-1) z = atan(y/y') -> dz/dx = cos(z)^2 - (n+0.5-x^2/4)*sin(z)^2
% 2-2-1-2) Direction fild for n = 1 over the full required range x in [0,10]
[X, Z] = meshgrid(linspace(0, 10, 30), linspace(-pi/2, pi/2, 25));
n_field = 1;
DZ = cos(Z).^2 - (n_field + 0.5 - X.^2/4).*sin(Z).^2;
DX = ones(size(DZ));
Lm = sqrt(DX.^2 + DZ.^2);
DX = DX./Lm; DZ = DZ./Lm;
figure('Name', 'Part 2: Direction Field (n=1)', 'NumberTitle', 'off');
quiver(X, Z, DX, DZ, 0.5);
grid on;
xlabel('x');
ylabel('z = arctan(y/y'')');
title('Direction Field for the z-Equation (n=1)');
axis([0 10 -pi/2 pi/2]);
% Turning point discussion (coefficient n+0.5-x^2/4 changes sign at x_t)
fprintf('2-2-1-2) Turning point x_t = 2*sqrt(n+0.5) where the coefficient (n+0.5-x^2/4) changes sign:\n');
for i = 1:length(n_vals)
    n = n_vals(i);
    x_turn = 2*sqrt(n+0.5);
    fprintf('    n = %2d  ->  x_t = %.4f\n', n, x_turn);
end
fprintf(['    Prediction: for x < x_t the coefficient is positive and y'''' = -(...)*y\n' ...
         '    behaves like a harmonic oscillator, so y(x) is OSCILLATORY near x = 0.\n' ...
         '    For x > x_t the coefficient becomes negative and y'''' = +(...)*y behaves\n' ...
         '    like an exponential term, so oscillations stop and the solution grows or\n' ...
         '    decays MONOTONICALLY beyond that point. x_t is exactly this critical\n' ...
         '    value where the qualitative behaviour changes.\n\n']);
% 2-2-2) Numerical solutions for n = 1 with both IC sets
% 3-2-2) Repeat for n = 5, 15 as well
figure('Name', 'Part 2: Parabolic Cylinder Solutions', 'NumberTitle', 'off');
colors = {'b', 'r', 'g'};
x_c2_all = cell(1, length(n_vals));
Y_c2_all = cell(1, length(n_vals));
for i = 1:length(n_vals)
    n = n_vals(i);
    ode_p2 = @(x, Y) [Y(2); -(n + 0.5 - x^2/4)*Y(1)];
    [x_c1, Y_c1] = ode45(ode_p2, xspan2, [1; 0]);   % y(0)=1, y'(0)=0
    [x_c2, Y_c2] = ode45(ode_p2, xspan2, [0; 1]);   % y(0)=0, y'(0)=1
    x_c2_all{i} = x_c2;
    Y_c2_all{i} = Y_c2;
    subplot(3,1,i);
    plot(x_c1, Y_c1(:,1), [colors{i} '-'], 'LineWidth', 2); hold on;
    plot(x_c2, Y_c2(:,1), [colors{i} '--'], 'LineWidth', 2);
    grid on;
    xlabel('x');
    ylabel('y(x)');
    title(sprintf('Parabolic Cylinder Solution for n = %d', n));
    legend('y(0)=1, y''(0)=0', 'y(0)=0, y''(0)=1');
end
% 3-2-2-c) Overlay comparison of the y(0)=0 solutions for n = 1, 5, 15
figure('Name', 'Part 2: Comparison of y(0)=0 Solutins Across n', 'NumberTitle', 'off');
hold on;
for i = 1:length(n_vals)
    plot(x_c2_all{i}, Y_c2_all{i}(:,1), colors{i}, 'LineWidth', 2, ...
         'DisplayName', sprintf('n = %d', n_vals(i)));
end
grid on;
xlabel('x');
ylabel('y(x)');
title('Comparison of y(0)=0, y''(0)=1 Solutions for n = 1, 5, 15');
legend('show');
fprintf(['3-2-2) Comparison and Analysis:\n' ...
         '  a) n = 5 vs n = 1: the turning point moves out to x_t = %.3f (from %.3f),\n' ...
         '     so the oscillatory region extends furthrr in x; within that region the\n' ...
         '     local "frequency" sqrt(n+0.5-x^2/4) is larger for n = 5, so its curve\n' ...
         '     oscillates faster (shorter local wavelength) than the n = 1 curve.\n' ...
         '  b) As n increases (1 -> 5 -> 15), x_t = 2*sqrt(n+0.5) grows, the\n' ...
         '     oscillatory region widens, the oscillations speed up, and the amplitude\n' ...
         '     grows more sharply once x exceeds x_t.\n' ...
         '  c) Comparing all three y(0) = 0 solutions:\n' ...
         '     Similarities:\n' ...
         '       - All three curves stsrt at the identical initial state\n' ...
         '         (y(0) = 0, y''(0) = 1).\n' ...
         '       - All exhibit local oscillatory behavior fir x < x_t and transition\n' ...
         '         to exponential/monotonic growth for x > x_t.\n' ...
         '     Differences:\n' ...
         '       - As n increases (1 -> 5 -> 15), the turning point\n' ...
         '         x_t = 2*sqrt(n+0.5) moves further to the right\n' ...
         '         (x_t = 2.45 -> 4.69 -> 7.87).\n' ...
         '       - Higher n values lead to significantly higher oscillation\n' ...
         '         frequency in the range 0 <= x < x_t.\n\n'], ...
         2*sqrt(5.5), 2*sqrt(1.5));
fprintf('--- Part 2 Execution Completed Successfully ---\n');
