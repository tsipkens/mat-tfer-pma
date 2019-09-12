
% MAIN_CHARGE  Script to investigate the effect of multiple chargining.
% Author:      Timothy Sipkens, 2019-06-25
%=========================================================================%


%-- Initialize script ----------------------------------------------------%
clear;
close all;

Rm = 10; % equivalent resolution of transfer functions (Reavell et al.)

m_star = 0.01e-18; % mass in kg (1 fg = 1e-18 kg)
m = linspace(1e-10,5,801).*m_star; % vector of mass

z_max = 4;
z_vec = [1,4];%1:z_max;
for zz=1:length(z_vec)
    z = z_vec(zz); % integer charge state
    disp(['Processing ',num2str(zz),' of ',num2str(length(z_vec)),'...']);
    
    rho_eff = 900; % effective density
    d = (6.*m./(rho_eff.*pi)).^(1/3);
        % specify mobility diameter vector with constant effective density
    
    prop = tfer_PMA.prop_PMA('Olfert'); % get properties of the CPMA
    % prop.omega_hat = 1; % NOTE: Uncomment for APM condition
    
    
    %=========================================================================%
    %-- Finite difference solution -------------------------------------------%
    tic;
    [tfer_FD(:,zz),sp,n{zz}] = tfer_PMA.tfer_FD(m_star,...
        m,d,z,prop,'Rm',Rm);
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
    %-- Diffusive transfer functions -----------------------------------------%
    %-- Method 1C --------------------------------%
    tic;
    tfer_1C_diff(:,zz) = ...
        tfer_PMA.tfer_1C_diff(m_star,m,d,z,prop,'Rm',Rm);
    t(12) = toc;
end

disp('Complete.');
disp(' ');


%=========================================================================%
%-- Plot different transfer functions with respect to m/m* ---------------%
m_plot = m./m_star;

figure(2);
plot(m_plot,tfer_1C_diff);
hold on;
plot(m_plot,min(tfer_FD,1),'k');
hold off;

ylim([0,1.2]);
% xlim(1+[-1.5/Rm,1.5/Rm]);

xlabel('m/m*')
ylabel('{\Lambda}')


