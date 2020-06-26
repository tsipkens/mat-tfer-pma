
% MAIN      Script used in plotting the cases from Olfert and Collings (2005).
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

prop = prop_pma('olfert-collings'); % get properties of the CPMA
prop.D = @(B) 1e-10.*ones(size(B));
prop.rho0 = rho_eff*pi/6; % copy mass-mobility relation info (only used to find Rm)
prop.Dm = 3;

prop_cpma = prop;
sp_cpma = get_setpoint(prop_cpma,'V',V,'omega',omega);

prop.omega_hat = 1;
sp = get_setpoint(prop,'V',V,'omega',omega);
    % get setpoint parameters



%=========================================================================%
%-- Finite difference solution -------------------------------------------%
[tfer_FD_w1] = tfer_FD(sp,m,d,1,prop); % for APM
[tfer_FD,sp_cpma] = ...
    tfer_FD(sp_cpma,m,d,1,prop_cpma); % for CPMA



%=========================================================================%
%-- Particle tracking approaches -----------------------------------------%
%-- Plug flow ------------------------------------------------------------%
%-- Method 1S ------------------------------%
[tfer_1S_w1] = tfer_1S_pb(sp,m,d,z,prop);



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


