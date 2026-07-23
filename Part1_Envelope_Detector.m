%% Differential Equations Final Project -- PART 1: Envelope Detector Circuit
% shima samadi , mahyas javaheri
% Term: 1404_05 Second Semester | Instructor: Dr. Shameli

clc;
clear;
close all;
fprintf('\n============================================================\n');
fprintf('PART 1: Envelope Detector Circuit\n');
fprintf('============================================================\n');
% Parametrs
C = 1e-6;    % 1 uF
R = 10e3;    % 10 kOhm
k = 0.05;    % Leakage factor (A/V^2)
V0 = 5;      % Initial voltage (V)
tspan1 = [0, 0.05];
% 1-2-1) Analyrical Solution
fprintf(['1-2-1) The ODE  dV/dt + (1/RC)V + (k/C)V^2 = 0  is a BERNOULLI\n' ...
         '  equation (y'' + P(t)y = Q(t)y^n, n = 2). Substituting u = V^(-1)\n' ...
         '  gives u'' - (1/RC)u = k/C, a first-order LINEAR ODE. Solving and\n' ...
         '  back-substituting (u = 1/V) yields the exact closed form:\n' ...
         '  V(t) = 1 / ( (k*R + 1/V0)*exp(t/(R*C)) - k*R )\n\n']);
% 1-2-2) Numerical Soltion using ode45
ode_p1 = @(t, V) -1/(R*C)*V - (k/C)*V^2;
[t1, V_num] = ode45(ode_p1, tspan1, V0);
% Analytical curve +ideal(k=0) refrence curve
V_anal  = 1 ./ ( (k*R + 1/V0) .* exp(t1./(R*C)) - k*R );
V_ideal = V0 * exp(-t1./(R*C));
% 1-2-3-1/2) Comparisn plot+ error plot
figure('Name', 'Part 1: Envelope Detector', 'NumberTitle', 'off');
subplot(2,1,1);
plot(t1, V_anal, 'r-', 'LineWidth', 2); hold on;
plot(t1, V_num, 'b--', 'LineWidth', 2);
plot(t1, V_ideal, 'g:', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Voltage V(t) [V]');
title('Envelope Detector Output Comparison');
legend('Analytical (k=0.05)', 'Numerical (ode45)', 'Ideal Case (k=0)');
subplot(2,1,2);
error1 = abs(V_anal - V_num);
semilogy(t1, error1, 'm-', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Absolute Error |Anal - Num|');
title('Numerical vs Analytical Error');
% 1-2-3-3) Discusion: efect of the nonlinear leakage term on discharge speed
thresh = V0/2;
t_half_real  = R*C * log( (1/thresh + k*R) / (k*R + 1/V0) );
t_half_ideal = R*C * log(2);
fprintf('1-2-3-3) Time for V(t) to fall to half its initial value (V0/2 = %.2f V):\n', thresh);
fprintf('    Real capacitor (k = 0.05) : t_1/2 = %.6e s\n', t_half_real);
fprintf('    Ideal capacitor (k = 0)   : t_1/2 = %.6e s\n', t_half_ideal);
fprintf(['    Because t_1/2(real) << t_1/2(ideal), thr nonlinear leakage term k*V^2\n' ...
         '    dominates the linear term V/R at these voltage levels and makes the\n' ...
         '    capacitor discharge MUCH faster thsn the ideal RC case. Physically,\n' ...
         '    "diagonal clipping distortion" happens when the RC discharge is too\n' ...
         '    slow to track a fast-dropping AM envelope; the extra nonlinear leakage\n' ...
         '    path speeds up the discharge (especially while V is large) and so tends\n' ...
         '    to REDUCE the risk of diagonal-clipping distortion compared to a purely\n' ...
         '    linear (ideal) detector with the ssme R and C.\n\n']);
fprintf('--- Part 1 Execution Completed Successfully ---\n');
