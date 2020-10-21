
% MAIN      Script used in plotting the cases from Ehara (1996).
% Author:   Timothy Sipkens, 2019-06-25
%=========================================================================%


%-- Initialize script ----------------------------------------------------%
clear;
close all;

V = 1000; % voltage to replicate Ehara et al.
omega = 1000*0.1047; % angular speed, converted from rpm to rad/s

e = 1.60218e-19; % electron charge [C]
m = linspace(260,340,601).*e; % vector of mass
% m = linspace(1e-10,6,601)*e; % vector of mass

z = 1; % integer charge state

rho_eff = 900; % effective density
d = (6.*m./(rho_eff.*pi)).^(1/3);
    % specify mobility diameter vector with constant effective density

prop = prop_pma('ehara'); % get properties of the CPMA
prop.m0 = rho_eff * pi / 6 * 1e-27; % copy mass-mobility relation info (only used to find Rm)
prop.Dm = 3;

sp = get_setpoint(prop,'V',V,'omega',omega);
    % get setpoint parameters

    
%=========================================================================%
%-- Finite difference solution -------------------------------------------%
tic;
[k_FD,n] = tfer_FD(sp,...
    m,d,1,prop);
t(1) = toc;


%=========================================================================%
%-- Transfer functions for different cases -------------------------------%
%-- Setup for centriputal force ------------------------------------------%
B = dm2zp(d,z,prop.T,prop.p);
tau = B.*m;
D = prop.D(B);
sig = sqrt(2.*prop.L.*D./prop.v_bar);
D0 = D.*prop.L/(prop.del^2*prop.v_bar); % dimensionless diffusion coeff.


%-- Particle tracking approaches -----------------------------------------%
%-- Plug flow ------------------------------------------------------------%
%-- Method 1S ------------------------------%
tic;
[k_1S,G0_1S] = tfer_1S(sp,m,d,z,prop);
t(1) = toc;

%-- Method 1S, Ehara et al. ----------------%
k_ehara = tfer_ehara(sp,m,d,z,prop);



%-- Parabolic flow -------------------------------------------------------%
%-- Method 1S ------------------------------%
tic;
[k_1S_pb,G0_1S_pb] = tfer_1S_pb(sp,m,d,z,prop);
t(2) = toc;



%=========================================================================%
%-- Plot different transfer functions with respect to m/m* ---------------%
m_plot = m./e;

figure(2);
plot(m_plot,k_1S);
hold on;
plot(m_plot,k_ehara);
plot(m_plot,k_1S_pb);
plot(m_plot,min(k_FD,1),'k');
hold off;

% ylim([0,1.2]);

xlabel('s')
ylabel('{\Lambda}')


