
% MAIN_VEC  A function to test vector evaluation of transfer functions.
%  AUTHOR: Timothy Sipkens, 2021-10-25

%-- Initialize script ----------------------------------------------------%
clear;
close all;

Rm = 3; % equivalent resolution of transfer functions (Reavell et al.)

m_star = logspace(-2, 0, 7) .* 1e-18; % mass in kg (1 fg = 1e-18 kg)
m = logspace(-2.3, 0.3, 801) .* 1e-18; % vector of mass

z = 1; % integer charge state

rho_eff100 = 900; % effective density
Dm = 3;
m100 = rho_eff100 * (pi * (100e-9)^3 / 6);
m0 = m100 * (1/100) ^ Dm;
d = 1e-9 .* (m ./ m0) .^ (1/Dm);
    % specify mobility diameter vector with constant effective density

prop = prop_pma('olfert'); % get properties of the CPMA
prop.omega_hat = 1;
prop.m0 = m0; % copy mass-mobility relation info (only used to find Rm)
prop.Dm = Dm;

% prop.omega_hat = 1; % NOTE: Uncomment for APM condition

sp = get_setpoint(prop, 'm_star', m_star, 'Rm', Rm);
    % get setpoint parameters


K0 = tfer_1C(sp, m, d, z, prop);
K1 = tfer_1C_diff(sp, m, d, z, prop);

%{
K2 = [];
for ii=1:length(sp)
    K2(ii,:) = tfer_FD(sp(ii), m, d, z, prop);
end
%}

figure(1);
semilogx(m, K0);
hold on;
semilogx(m, K1, 'k--');
% semilogx(m, K2, 'r--');
hold off;


