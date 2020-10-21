
% MAIN    Script used in plotting different transfer functions.
%         This is used to evaluate all of the available methods of
%         determining the transfer function and replicates the plots
%         in the associated journal article (details of that work are
%         included in the README).
%         associated
% Author: Timothy Sipkens, 2019-06-25
%=========================================================================%


%-- Initialize script ----------------------------------------------------%
clear;
close all;

Rm = 10; % equivalent resolution of transfer functions (Reavell et al.)

m_star = 0.01e-18; % mass in kg (1 fg = 1e-18 kg)
m = linspace(0.8,1.2,601).*m_star; % vector of mass

z = 1; % integer charge state

rho_eff = 900; % effective density
Dm = 3;
m0 = pi * rho_eff / 6 * 100 ^ (3-Dm);
d = (6.*m./(rho_eff.*pi)).^(1/Dm);
    % specify mobility diameter vector with constant effective density

prop = prop_pma('olfert'); % get properties of the CPMA
prop.m0 = rho_eff*pi/6; % copy mass-mobility relation info (only used to find Rm)
prop.Dm = Dm;

% prop.omega_hat = 1; % NOTE: Uncomment for APM condition

sp = get_setpoint(prop,'m_star',m_star,'Rm',Rm);
    % get setpoint parameters


%=========================================================================%
%-- Finite difference solution -------------------------------------------%
tic;
[k_FD,n] = tfer_FD(sp,...
    m,d,z,prop);
t(1) = toc;


%=========================================================================%
%-- Particle tracking approaches -----------------------------------------%
%-- Plug flow ------------------------------------------------------------%
%-- Method 1S ------------------------------%
tic;
[k_1S,G0_1S] = tfer_1S(sp,m,d,z,prop);
t(2) = toc;

%-- Method 1S, Ehara et al. ----------------%
k_ehara = tfer_ehara(sp,m,d,z,prop);

%-- Method 1C ------------------------------%
tic;
[k_1C,G0_1C] = tfer_1C(sp,m,d,z,prop);
t(3) = toc;

%-- Method 2S ------------------------------%
tic;
[k_2S,G0_2S] = tfer_2S(sp,m,d,z,prop);
t(4) = toc;

%-- Method 2C ------------------------------%
tic;
[k_2C,G0_2C] = tfer_2C(sp,m,d,z,prop);
t(5) = toc;

%-- Method W1 ------------------------------%
if prop.omega_hat==1
    tic;
    [k_W1,G0_W1] = tfer_W1(sp,m,d,z,prop);
    t(6) = toc;
end

%-- Method GE ------------------------------%
tic;
[k_GE,G0_GE] = tfer_GE(sp,m,d,z,prop);
t(7) = toc;


%-- Parabolic flow -------------------------------------------------------%
%-- Method 1S ------------------------------%
tic;
[k_1S_pb,G0_1S_pb] = tfer_1S_pb(sp,m,d,z,prop);
t(8) = toc;

%-- Method 1C ------------------------------%
tic;
[k_1C_pb,G0_1C_pb] = tfer_1C_pb(sp,m,d,z,prop);
t(9) = toc;

%-- Method W1 ------------------------------%
if prop.omega_hat==1
    tic;
    [k_W1_pb,G0_W1_pb] = tfer_W1_pb(sp,m,d,z,prop);
    t(10) = toc;
end


%-- Diffusive transfer functions -----------------------------------------%
%-- Method 1S --------------------------------%
tic;
k_1S_diff = tfer_1S_diff(sp,m,d,z,prop);
t(11) = toc;

%-- Method 1C --------------------------------%
tic;
k_1C_diff = tfer_1C_diff(sp,m,d,z,prop);
t(12) = toc;

%-- Method 2S --------------------------------%
tic;
k_2S_diff = tfer_2S_diff(sp,m,d,z,prop);
t(13) = toc;

%-- Method 2C --------------------------------%
tic;
k_2C_diff = tfer_2C_diff(sp,m,d,z,prop);
t(14) = toc;

%-- Method W1 --------------------------------%
if prop.omega_hat==1
    tic;
    k_W1_diff = tfer_W1_diff(sp,m,d,z,prop);
    t(15) = toc;
end

%-- Method GE --------------------------------%
tic;
k_GE_diff = tfer_GE_diff(sp,m,d,z,prop);
t(16) = toc;


%-- Triangle approx. -----------------------%
tic;
k_tri = tfer_tri(sp,m,d,z,prop);
t(18) = toc;



%=========================================================================%
%-- Plot different transfer functions with respect to m/m* ---------------%
m_plot = m./m_star;

figure(2);
% plot(m_plot,k_1S);
% hold on;
% plot(m_plot,k_ehara);
% plot(m_plot,k_1S_diff);
% hold on;
% plot(m_plot,k_1S_pb);
plot(m_plot,k_1C);
hold on;
plot(m_plot,k_1C_diff);
plot(m_plot,k_1C_pb);
% plot(m_plot,k_2S);
% plot(m_plot,k_2S_diff);
% plot(m_plot,k_2C);
% plot(m_plot,k_2C_diff);
% plot(m_plot,k_W1,'r');
% plot(m_plot,k_W1_diff,'r');
% plot(m_plot,k_W1_pb);
% plot(m_plot,k_GE);
% plot(m_plot,k_tri);
plot(m_plot,min(k_FD,1),'k');
hold off;

ylim([0,1.2]);
xlim(1+[-1.5/Rm,1.5/Rm]);

xlabel('m/m*')
ylabel('{\Lambda}')



