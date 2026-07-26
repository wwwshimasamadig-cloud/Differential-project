%% Differential Equations Final Project -- PART 5: Circuit Analysis with Matrices
% shima samadi , mahyas javaheri
% Term: 1404_05 Second Semester | Instructor: Dr. Shameli

clc;
clear;
close all;
fprintf('\n============================================================\n');
fprintf('PART 5: Matrix Analysis of RLC Circuit\n');
fprintf('============================================================\n');
L5 = 1; C5 = 0.4; R5 = 3;
tspan5 = [0, 5];
X0 = [5; 0];
% 5-2-1) Matrix Formoulation
A = [ 0,      1/C5;
     -1/L5, -R5/L5 ];
fprintf('State Space Matrix A:\n');
disp(A);
[V_eig, D_eig] = eig(A);
eigenvalues = diag(D_eig);
fprintf('Eigenvalues lambda_1, lambda_2:\n');
disp(eigenvalues);
% 5-2-1-2) Analytical solution from eigenvalues -1.5 +/- 0.5i and X0 = [5;0]
% VC(t) = exp(-1.5*t)*(5*cos(0.5*t) + 15*sin(0.5*t))
% i(t)  = C5*VC'(t) = -10*exp(-1.5*t)*sin(0.5*t)
t_anal5 = linspace(tspan5(1), tspan5(2), 200);
VC_analytical = exp(-1.5*t_anal5) .* (5*cos(0.5*t_anal5) + 15*sin(0.5*t_anal5));
i_analytical  = -10*exp(-1.5*t_anal5) .* sin(0.5*t_anal5);
% 5-2-2) Numerical Solver Integration via ode45
matrix_ode = @(t, X) A*X;
[t_num5, X_num5] = ode45(matrix_ode, t_anal5, X0);
% 5-2-3)Superimposed Graphics Comparison
figure('Name', 'Part 5: Matrix RLC Transient Response', 'NumberTitle', 'off');
subplot(2,1,1);
plot(t_anal5, VC_analytical, 'r-', 'LineWidth', 2.5); hold on;
plot(t_num5, X_num5(:,1), 'b--', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Capacitor Voltage V_C(t) [V]');
title('Capacitor Voltage Comparison');
legend('Analytical Solution', 'Numerical Solotion (ode45)');
subplot(2,1,2);
plot(t_anal5, i_analytical, 'r-', 'LineWidth', 2.5); hold on;
plot(t_num5, X_num5(:,2), 'b--', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Circuit Current i(t) [A]');
title('Inductor Current Comparison');
legend('Analytical Solution', 'Numerical Solution (ode45)');
% 3-2-5-2) Analysis: matching accuracy +damping classification
err_VC = max(abs(VC_analytical - X_num5(:,1)'));
err_i  = max(abs(i_analytical - X_num5(:,2)'));
fprintf('3-2-5-2) Max |analytical - numerical| error: VC -> %.3e V , i -> %.3e A\n', err_VC, err_i);
fprintf(['    -> The hand (analytical) solution matches the ode45 simulation to\n' ...
         '    numerical precision.\n' ...
         '    The eigenvalues lambda = %.2f +/- %.2fi are COMPLEX with a NEGATIVE\n' ...
         '    real part, so the system is UNDERDAMPED: the response is a decaying\n' ...
         '    sinusoidal oscillation, with decay rate |Re(lambda)| = %.2f 1/s and\n' ...
         '    damped natural frequency |Im(lambda)| = %.2f rad/s.\n\n'], ...
         real(eigenvalues(1)), abs(imag(eigenvalues(1))), ...
         abs(real(eigenvalues(1))), abs(imag(eigenvalues(1))));
fprintf('--- Part 5 Execution Completed Successfully ---\n');