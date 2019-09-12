
% MAIN      Script used in plotting different transfer functions.
% Author:   Timothy Sipkens, 2019-06-25
%=========================================================================%


%-- Initialize script ----------------------------------------------------%
clear;
close all;

V = 100; % voltage to replicate Ehara et al.
omega = 5000*0.1047; % angular speed, converted from rpm to rad/s

e = 1.60218e-19; % electron charge [C]
m = linspace(0.4,0.7,601).*1e-18; % vector of mass
% m = linspace(1e-10,6,601)*e; % vector of mass

z = 1; % integer charge state

d = 100e-9.*ones(size(m));
    % specify mobility diameter vector with constant effective density

prop = tfer_PMA.prop_PMA('Kuwata'); % get properties of the CPMA
prop.D = @(B) 1e-10.*ones(size(B));



%=========================================================================%
%-- Finite difference solution -------------------------------------------%
[tfer_FD,sp] = tfer_PMA.tfer_FD([],...
    m,d,1,prop,'V',V,'omega',omega);



%=========================================================================%
%-- Transfer functions for different cases -------------------------------%
%-- Setup for centriputal force ------------------------------------------%
B = tfer_PMA.dm2zp(d,z,prop.T,prop.p);
tau = B.*m;

%-- Particle tracking approaches -----------------------------------------%
%-- Plug flow ------------------------------------------------------------%
%-- Method 1S ------------------------------%
[tfer_1S] = ...
    tfer_PMA.tfer_1S([],m,d,z,prop,'V',V,'omega',omega);
[tfer_1S_pb] = ...
    tfer_PMA.tfer_1S_pb([],m,d,z,prop,'V',V,'omega',omega);



%=========================================================================%
%-- Plot different transfer functions with respect to m/m* ---------------%
m_plot = m;

figure(2);
plot(m_plot,min(tfer_FD,1));
hold on;
plot(m_plot,min(tfer_1S,1));
plot(m_plot,min(tfer_1S_pb,1));
hold off;

% ylim([0,1.2]);

xlabel('s')
ylabel('{\Lambda}')


