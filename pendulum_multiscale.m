%% pendulum_multiscale.m
%
% Numerical solution of the pendulum equation
%
%       theta'' + sin(theta) = 0,   theta(0) = eps,  theta'(0) = 0
%
% compared against the leading-order multiple-scales approximation
% derived in Part 2 of the project:
%
%       theta(t) ~ eps * cos( t - (eps^2/16) t )
%                = eps * cos( (1 - eps^2/16) t )
%
% and the corresponding approximate period
%
%       T ~ 2*pi / (1 - eps^2/16)  ~  2*pi*(1 + eps^2/16)
%
% -------------------------------------------------------------------

clear; close all; clc;

eps_list = [0.3,0.6];
tmin = 0;
tmax = 80;                
Npts = 4000;              

opts = odeset('RelTol', 1e-10, 'AbsTol', 1e-12);
 
for eps = eps_list
 
    % numerical solution (sol is a struct; deval evaluates it at any t)
    sol = ode45(@(t,y) [y(2); -sin(y(1))], [0 80], [eps; 0], opts);
 
    t = linspace(tmin, tmax, Npts);
    theta_num  = deval(sol, t, 1);
    theta_asym = eps*cos((1 - eps^2/16)*t);
 
    figure; plot(t, theta_num, t, theta_asym, '--');
    xlabel('t'); ylabel('\theta'); legend('numerical', 'asymptotic');
    title(['\epsilon = ' num2str(eps)]);
    grid on;
 
    % period: theta first hits zero at T/4 (start is at rest at theta = eps)
    T_num  = 4*fzero(@(t) deval(sol, t, 1), [1 2]);
    T_asym = 2*pi/(1 - eps^2/16);
 
    fprintf('eps = %.2f:  T_num = %.6f,  T_asym = %.6f\n', eps, T_num, T_asym);
end
