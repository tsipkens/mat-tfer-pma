
% MAIN      Script used in plotting different transfer functions.
% Author:   Timothy Sipkens, 2019-06-25
%=========================================================================%


%-- Initialize script ----------------------------------------------------%
clear;
close all;

Rm = 10; % equivalent resolution of transfer functions (Reavell et al.)

m_star = 0.01e-18; % mass in kg (1 fg = 1e-18 kg)
m = linspace(0.8,1.2,601).*m_star; % vector of mass

z = 1; % integer charge state

rho_eff = 900; % effective density
d = (6.*m./(rho_eff.*pi)).^(1/3);
    % specify mobility diameter vector with constant effective density

prop = tfer_PMA.prop_PMA('Olfert'); % get properties of the CPMA
prop.omega_hat = 1; % NOTE: Uncomment for APM condition


%=========================================================================%
%-- Finite difference solution -------------------------------------------%
tic;
[tfer_FD,sp,n] = tfer_PMA.tfer_FD(m_star,...
    m,d,1,prop,'Rm',Rm);
t(1) = toc;


%=========================================================================%
%-- Transfer functions for different cases -------------------------------%
%-- Setup for centriputal force ------------------------------------------%
if ~exist('d','var')
    B = tfer_PMA.mp2zp(m,z,prop.T,prop.p);
else
    B = tfer_PMA.dm2zp(d,z,prop.T,prop.p);
end
tau = B.*m;
D = prop.D(B);
sig = sqrt(2.*prop.L.*D./prop.v_bar);
D0 = D.*prop.L/(prop.del^2*prop.v_bar); % dimensionless diffusion coeff.


%-- Particle tracking approaches -----------------------------------------%
%-- Plug flow ------------------------------------------------------------%
%-- Method 1S ------------------------------%
tic;
[tfer_1S,G0_1S] = tfer_PMA.tfer_1S(m_star,m,d,z,prop,'Rm',Rm);
t(2) = toc;

%-- Method 1S, Ehara et al. ----------------%
tfer_Ehara = tfer_PMA.tfer_Ehara(m_star,m,d,z,prop,'Rm',Rm);

%-- Method 1C ------------------------------%
tic;
[tfer_1C,G0_1C] = tfer_PMA.tfer_1C(m_star,m,d,z,prop,'Rm',Rm);
t(3) = toc;

%-- Method 2S ------------------------------%
tic;
[tfer_2S,G0_2S] = tfer_PMA.tfer_2S(m_star,m,d,z,prop,'Rm',Rm);
t(4) = toc;

%-- Method 2C ------------------------------%
tic;
[tfer_2C,G0_2C] = tfer_PMA.tfer_2C(m_star,m,d,z,prop,'Rm',Rm);
t(5) = toc;

%-- Method W1 ------------------------------%
if prop.omega_hat==1
    tic;
    [tfer_W1,G0_W1] = tfer_PMA.tfer_W1(m_star,m,d,z,prop,'Rm',Rm);
    t(6) = toc;
end

%-- Method GE ------------------------------%
tic;
[tfer_GE,G0_GE] = tfer_PMA.tfer_GE(m_star,m,d,z,prop,'Rm',Rm);
t(7) = toc;


%-- Parabolic flow -------------------------------------------------------%
%-- Method 1S ------------------------------%
tic;
[tfer_1S_pb,G0_1S_pb] = tfer_PMA.tfer_1S_pb(m_star,m,d,z,prop,'Rm',Rm);
t(8) = toc;

%-- Method 1C ------------------------------%
tic;
[tfer_1C_pb,G0_1C_pb] = tfer_PMA.tfer_1C_pb(m_star,m,d,z,prop,'Rm',Rm);
t(9) = toc;

%-- Method W1 ------------------------------%
if prop.omega_hat==1
    tic;
    [tfer_W1_pb,G0_W1_pb] = tfer_PMA.tfer_W1_pb(m_star,m,d,z,prop,'Rm',Rm);
    t(10) = toc;
end


%-- Diffusive transfer functions -----------------------------------------%
%-- Method 1S --------------------------------%
tic;
tfer_1S_diff = tfer_PMA.tfer_1S_diff(m_star,m,d,z,prop,'Rm',Rm);
t(11) = toc;

%-- Method 1C --------------------------------%
tic;
tfer_1C_diff = tfer_PMA.tfer_1C_diff(m_star,m,d,z,prop,'Rm',Rm);
t(12) = toc;

%-- Method 2S --------------------------------%
tic;
tfer_2S_diff = tfer_PMA.tfer_2S_diff(m_star,m,d,z,prop,'Rm',Rm);
t(13) = toc;

%-- Method 2C --------------------------------%
tic;
tfer_2C_diff = tfer_PMA.tfer_2C_diff(m_star,m,d,z,prop,'Rm',Rm);
t(14) = toc;

%-- Method W1 --------------------------------%
if prop.omega_hat==1
    tic;
    tfer_W1_diff = tfer_PMA.tfer_W1_diff(m_star,m,d,z,prop,'Rm',Rm);
    t(15) = toc;
end

%-- Method GE --------------------------------%
tic;
tfer_GE_diff = tfer_PMA.tfer_GE_diff(m_star,m,d,z,prop,'Rm',Rm);
t(16) = toc;


%-- Triangle approx. -----------------------%
tic;
tfer_tri = tfer_PMA.tfer_tri(m_star,m,d,z,prop,'Rm',Rm);
t(18) = toc;



%=========================================================================%
%-- Plot different transfer functions with respect to m/m* ---------------%
m_plot = m./m_star;

figure(2);
% plot(m_plot,tfer_1S);
% hold on;
% plot(m_plot,tfer_Ehara);
% plot(m_plot,tfer_1S_diff);
% hold on;
% plot(m_plot,tfer_1S_pb);
plot(m_plot,tfer_1C);
hold on;
plot(m_plot,tfer_1C_diff);
plot(m_plot,tfer_1C_pb);
% plot(m_plot,tfer_2S);
% plot(m_plot,tfer_2S_diff);
% plot(m_plot,tfer_2C);
% plot(m_plot,tfer_2C_diff);
% plot(m_plot,tfer_W1,'r');
% plot(m_plot,tfer_W1_diff,'r');
% plot(m_plot,tfer_W1_pb);
% plot(m_plot,tfer_GE);
% plot(m_plot,tfer_tri);
plot(m_plot,min(tfer_FD,1),'k');
hold off;

ylim([0,1.2]);
xlim(1+[-1.5/Rm,1.5/Rm]);

xlabel('m/m*')
ylabel('{\Lambda}')

%{
%=========================================================================%
%-- Bar plot of error ----------------------------------------------------%
chi_sq = [];
mean_sq_err = [];

vec = {'1S','1C','2S','2C','W1','GE','1S_pb','1C_pb','W1_pb',...
    '1S_diff','1C_diff','2S_diff','2C_diff','W1_diff','GE_diff','tri'};
for ii=1:length(vec)
    if ~and(prop.omega_hat~=1,~isempty(strfind(vec{ii},'W')))
        ind_nz = or(eval(['tfer_',vec{ii},'>1e-5']),tfer_FD>1e-5);
        chi_sq{ii} = (eval(['tfer_',vec{ii},'(ind_nz)'])-tfer_FD(ind_nz)).^2;
        mean_sq_err(ii) = mean(chi_sq{ii});
    end
end

figure(3);
semilogy(mean_sq_err,'o-');
ylim([1e-6,1e-2]);
%}

