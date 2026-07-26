%% Differential Equations Final Project -- PART 3: Time-Varying Conductivity RLC Circuit
% shima samadi , mahyas javaheri
% Term: 1404_05 Second Semester | Instructor: Dr.Shameli

clc;
clear;
close all;
fprintf('\n============================================================\n');
fprintf('PART 3: Time-Varying Conductivity RLC Circuit\n');
fprintf('============================================================\n');
L_val = 1;        % 1 H
C_val = 1e-3;     % 1 mF
G0 = 0.02;        % Initial conductance (S)
alpha = 0.5;      % Rate of conductance increase (S/s)
tspan3 = [0, 0.005];
v0_cap = 5;
i0_ind = 0;
% Reference ("ground truth " ) numerical solution
t_series = linspace(tspan3(1), tspan3(2), 400);
ode_p3 = @(t, Y) [ -((G0 + alpha*t)*Y(1) + Y(2))/C_val; Y(1)/L_val ];
[t_num3, Y_num3] = ode45(ode_p3, t_series, [v0_cap; i0_ind]);
% 2-2-3) Reusable iterative coefficient generator: power_series_coeffs()
% (function definition at the end of this file) takes (L,C,G0,alpha,v0,i0,N)
% and returns the power-series coefficients C_0 ... C_N for v(t).
N_terms = 30;
C_coeff = power_series_coeffs(L_val, C_val, G0, alpha, v0_cap, i0_ind, N_terms);
v_series = zeros(size(t_series));
for n = 0:N_terms
    v_series = v_series + C_coeff(n+1) * (t_series.^n);
end
figure('Name', 'Part 3: Time-Varying Sensor Circuit', 'NumberTitle', 'off');
plot(t_series, v_series, 'r-', 'LineWidth', 2); hold on;
plot(t_num3, Y_num3(:,1), 'b--', 'LineWidth', 2);
grid on;
xlabel('Time (s)');
ylabel('Capacitor Voltage v(t) [V]');
title('Power Series Solution vs. ode45 Reference');
legend(sprintf('Power Series (N=%d)', N_terms), 'ode45 Solution');
% 3-2-3-1) Convergence radius: repeat for several N and inspect the eror
N_list = [5, 10, 20, 30, 50];
figure('Name', 'Part 3: Convergence of the Power Series', 'NumberTitle', 'off');
hold on;
for j = 1:length(N_list)
    Nj = N_list(j);
    Cj = power_series_coeffs(L_val, C_val, G0, alpha, v0_cap, i0_ind, Nj);
    vj = zeros(size(t_series));
    for n = 0:Nj
        vj = vj + Cj(n+1) * (t_series.^n);
    end
    err_j = abs(vj - Y_num3(:,1)');
    semilogy(t_series, err_j, 'LineWidth', 1.5, 'DisplayName', sprintf('N = %d', Nj));
end
set(gca, 'YScale', 'log');
grid on;
xlabel('Time (s)');
ylabel('|Power series - ode45|');
title('Truncation Error vs. Time for Different Series Lengths N');
legend('show');
fprintf(['3-2-3-1) As the convergence plot shows, larger N keeps the power series\n' ...
         '  accurate over a longer stretch of the [0, 5 ms] window before the\n' ...
         '  truncation error turns sharply upward; that turning point marks the\n' ...
         '  effective radius of convergence for that N. More terms push the\n' ...
         '  breakdown point later, but for anu fixed N the series eventually\n' ...
         '  diverges once t grows past that N-dependent radius.\n\n']);
fprintf(['3-2-3-2) Physical discussion: the sensor conductance G(t) = G0 + alpha*t\n' ...
         '  increases linearly with time, which increases the damping term (G(t)/C)\n' ...
         '  in the governing ODE. So the effective damping of the circuit rises over\n' ...
         '  time: any oscillatory tendency (the roots ar t=0 are complex/underdamped)\n' ...
         '  gets suppressed increasingly fast, and the circuit behaves more and more\n' ...
         '  like the capacitor bwing progressively short-circuited by the growing\n' ...
         '  sensoor conductance, driving v(t) toward zero faster than a constant-G\n' ...
         '  RLC circuit would.\n\n']);
fprintf('--- Part 3 Execution Completed Successfully ---\n');
%========================================================================
% LOCAL FUNCTINS
%========================================================================
function C_coeff = power_series_coeffs(L_val, C_val, G0, alpha, v0, i0, N_terms)
    % power_series_coeffs: iterative generator for the Taylor/power-series
    % coefficients of v(t) = sum_{n=0}^{N_terms} C_n * t^n, where v(t)
    % satisfies:
    %   v'' + ((G0 + alpha*t)/C) * v' + ((L*alpha + 1)/(L*C)) * v = 0
    % (the time-varing-conductance RLC circuit of Part 3), with initial
    % conditions v(0) = v0 and C*v'(0) +G0*v0 + i0 = 0 (KCL at t = 0).
    % Inputs : L_val, C_val, G0, alpha- circuit parameters
    %          v0, i0                   - initial capacitor voltage / inductor current
    %          N_terms                 - highest power of t to compute
    % Output : C_coeff (1 x (N_terms+1)) - coefficients C_0 ... C_{N_terms}
    C_coeff = zeros(1, N_terms + 1);
    C_coeff(1) = v0;                              % C0 = v(0)
    C_coeff(2) = -(G0*v0 + i0) / C_val;         % C1 = v'(0)
    for n = 0:(N_terms - 2)
        term1 = (G0 / C_val) * (n + 1) * C_coeff(n + 2);
        term2 = (alpha / C_val) * n * C_coeff(n + 1);
        term_v = ((L_val*alpha + 1) / (L_val*C_val)) * C_coeff(n + 1);
        C_coeff(n + 3) = -(term1 + term2 + term_v) / ((n + 2) * (n + 1));
    end
end
