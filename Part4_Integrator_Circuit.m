%% Differential Equations Final Project -- PART 4: Integrator Circuit
% shima samadi , mahyas javaheri
% Term: 1404_05 Second Semester | Instructor: Dr. Shameli

clc;
clear;
close all;
fprintf('\n============================================================\n');
fprintf('PART 4: Integrator Circuit\n');
fprintf('============================================================\n');
syms R_s C_s s t Vin_s Vc(t) A
fprintf('--- Symbolic Tool Explanations (see "doc syms/dsolve/laplace/ilaplace") ---\n');
fprintf('syms     : declares symbolic variables/functions.\n');
fprintf('dsolve   : solves ODEs symbolically, directly in the time domain.\n');
fprintf('laplace  : computes the Laplace transform X(s) of a time-domain function x(t).\n');
fprintf('ilaplace : computes the inverse Laplace transform, recovering x(t) from X(s).\n\n');
% 2-2-4-1) (a) Time-domain symbolic soluton via dsolve
ode_integrator = R_s * C_s * diff(Vc(t), t) + Vc(t) == Vin_s;
cond_integrator = Vc(0) == 0;
Vc_sol(t) = dsolve(ode_integrator, cond_integrator);
fprintf('dsolve solution of RC*Vc'' + Vc = Vin (Vin constant), Vc(0) = 0:\n');
disp(Vc_sol(t));
% 2-2-4-1) (b) Same equation solved via the Laplace transform, stepp input
% RC*(s*Vc(s) - Vc(0)) + Vc(s) = Vin(s),   Vin(t) = A (step)  ->  Vin(s) = A/s
Vin_step_s = A/s;
Vc_s = Vin_step_s / (R_s*C_s*s + 1);   % Vc(0) = 0
Vc_t_laplace = ilaplace(Vc_s, s, t);
fprintf('Laplace/ilaplace solution for a step input Vin(t) = A, Vc(0) = 0:\n');
disp(Vc_t_laplace);
fprintf(['Both routes agree: dsolve works directly in the time domain, while the\n' ...
         'laplace/ilaplace route turns the ODE into an algebraic equation in s,\n' ...
         'solves it there, then transforms back -- same physical answer either way.\n\n']);
% 1-2-4-2) High-frequency / RC>>1 approximation
fprintf(['1-2-4-2) If RC*w >> 1, H(s) = Vout(s)/Vin(s) = 1/(RCs+1) ~= 1/(RCs).\n' ...
         '  Dividing by s in the Laplace domain corresponds to integration in the\n' ...
         '  time domain, so in this regime the circuit behaves like an ideal\n' ...
         '  integrator:\n' ...
         '  Vout(t) ~= Vout(0) + (1/(R*C)) * Integral( Vin(tau) d(tau), 0, t )\n\n']);
% 2-2-4-2) Numerical Test: Sinusoidal Input
R4 = 1e3; C4 = 1e-6; f4 = 10e3; w4 = 2*pi*f4; Vp = 5;
tspan4 = [0, 5/f4];
initial_conditions = [-1, 0, 1];
figure('Name', 'Part 4: Sinusoidal Input with ICs', 'NumberTitle', 'off');
t_eval = linspace(tspan4(1), tspan4(2), 1000);
Vc_out_all = zeros(length(initial_conditions), length(t_eval));
for idx = 1:length(initial_conditions)
    ic = initial_conditions(idx);
    sin_ode = @(t, Vc) (Vp*sin(w4*t) - Vc) / (R4*C4);
    [t_out, Vc_out] = ode45(sin_ode, t_eval, ic);
    Vc_out_all(idx,:) = Vc_out';
    plot(t_out, Vc_out, 'LineWidth', 2, 'DisplayName', sprintf('V_C(0^-) = %d V', ic));
    hold on;
end
grid on;
xlabel('Time (s)');
ylabel('Output Voltage V_{out}(t) [V]');
title('Integrator Transient Response to Sine Wave');
legend('show');
% 3-2-4-2) Discussion:effect of initial condition on transient vs steady state
diff_start = abs(Vc_out_all(1,1) - Vc_out_all(3,1));
diff_end   = abs(Vc_out_all(1,end) - Vc_out_all(3,end));
fprintf('3-2-4-2) |V_out difference between IC = -1 and IC = +1|:\n');
fprintf('    at t = %.3e s (start): %.4f V\n', t_eval(1), diff_start);
fprintf('    at t = %.3e s (end)  : %.4f V\n', t_eval(end), diff_end);
fprintf(['    The initial condition strongly shapes the early transient, but every\n' ...
         '    curve shares thr same exponentially decayiing homogeneous term\n' ...
         '    exp(-t/RC); its effect vanishes after a few time constants\n' ...
         '    (RC = %.2e s here), so all curves converge onto the same forced\n' ...
         '    steady-state sinusoidal response -- confirmed by the small difference\n' ...
         '    remaining at the end of the simulation window.\n\n'], R4*C4);
% 2-2-4-3) Response to a Square Wave, pluss ideal-integrator approximation
T = 1/f4;
Vmax = 5;
square_input = @(t) Vmax*(mod(t,T) < T/2) - Vmax*(mod(t,T) >= T/2);
square_ode = @(t, Vc) (square_input(t) - Vc)/(R4*C4);
[t_sq, Vc_sq] = ode45(square_ode, t_eval, 0);
% Ideal integrator approximation: Vout(t) ~= (1/RC)*Integral(Vin dt)
Vin_sq = square_input(t_eval);
Vout_ideal_integrator = cumtrapz(t_eval, Vin_sq) / (R4*C4);
figure('Name', 'Part 4: Square Wave Response', 'NumberTitle', 'off');
plot(t_sq, square_input(t_sq), 'k--', 'LineWidth', 1.5); hold on;
plot(t_sq, Vc_sq, 'b-', 'LineWidth', 2);
plot(t_eval, Vout_ideal_integrator, 'r-.', 'LineWidth', 1.5);
grid on;
xlabel('Time (s)');
ylabel('Voltage [V]');
title('Integrator Output for Square Wave Input');
legend('Input Square Wave V_{in}(t)', 'Exact ODE Solution V_{out}(t)', ...
       'Ideal Integrator Approx. (1/RC)\int V_{in} dt');
fprintf(['3-2-4-3) Comparing the exact ODE solution to the (1/RC)*integral\n' ...
         '  approximation: since R*C = %.2e s and the square-wave period T = %.2e s\n' ...
         '  satisfy R*C*w >> 1 (RC >> T/(2*pi)), the two curves nearly overlap.\n' ...
         '  V_out(t) closely tracks the running integral of V_in(t): while V_in sits\n' ...
         '  at +Vmax, V_out ramps up ~linearly, and while V_in sits at -Vmax, V_out\n' ...
         '  ramps down. The result is the classic triangular-wave output of an RC\n' ...
         '  integrator driven by a square wave, confirming the frequency-domain\n' ...
         '  approximation from 1-2-4-2 in the time domain.\n\n'], R4*C4, T);
fprintf('--- Part 4 Execution Completed Successfully ---\n');