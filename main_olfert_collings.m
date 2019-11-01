
% MAIN      Script used in plotting different transfer functions.
% Author:   Timothy Sipkens, 2019-06-25
%=========================================================================%


%-- Initialize script ----------------------------------------------------%
clear;
close all;

V = 100; % voltage to replicate Ehara et al.
omega = 2500*0.1047; % angular speed, converted from rpm to rad/s

e = 1.60218e-19; % electron charge [C]
m = linspace(4,6,601)*e; % vector of mass
% m = linspace(1e-10,6,601)*e; % vector of mass

z = 1; % integer charge state

rho_eff = 900; % effective density
d = (6.*m./(rho_eff.*pi)).^(1/3);
    % specify mobility diameter vector with constant effective density

prop = tfer_pma.prop_pma('Olfert-Collings'); % get properties of the CPMA
prop.D = @(B) 1e-10.*ones(size(B));
omega_hat = prop.omega_hat; % only valid for CPMA



%=========================================================================%
%-- Finite difference solution -------------------------------------------%
prop.omega_hat = 1;
[tfer_FD_w1,sp] = tfer_pma.tfer_FD([],...
    m,d,1,prop,'V',V,'omega',omega);

prop.omega_hat = omega_hat;
[tfer_FD,sp_cpma] = tfer_pma.tfer_FD([],...
    m,d,1,prop,'V',V,'omega',omega);



%=========================================================================%
%-- Transfer functions for different cases -------------------------------%
%-- Setup for centriputal force ------------------------------------------%
prop = tfer_pma.prop_pma('Olfert-Collings'); % get properties of the CPMA
B = tfer_pma.dm2zp(d,z,prop.T,prop.p);
tau = B.*m;

%-- Particle tracking approaches -----------------------------------------%
%-- Plug flow ------------------------------------------------------------%
%-- Method 1S ------------------------------%
prop.omega_hat = 1;
[tfer_1S_w1] = ...
    tfer_pma.tfer_1S_pb([],m,d,z,prop,'V',V,'omega',omega);



%=========================================================================%
%-- Plot different transfer functions with respect to m/m* ---------------%
m_plot = m./e;

figure(2);
plot(m_plot,min(tfer_FD,1));
hold on;
plot(m_plot,min(tfer_FD_w1,1));
plot(m_plot,min(tfer_1S_w1,1));
hold off;

% ylim([0,1.2]);

xlabel('s')
ylabel('{\Lambda}')


