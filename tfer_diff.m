
% TFER_DIFF     General function for diffusion, only requiring a different G0.
% Author:       Timothy Sipkens, 2018-12-27
% 
% Inputs:
%   sp          Structure defining various setpoint parameters 
%               (e.g. m_star, V). Use 'get_setpoint' method to generate 
%               this structure.
%   m           Particle mass
%   d           Particle mobility diameter
%   z           Integer charge state
%   prop        Device properties (e.g. classifier length)
%
% Outputs:
%   Lambda      Transfer function
%   G0          Function mapping final to initial radial position
%=========================================================================%

function [Lambda] = tfer_diff(G0, D, prop)

sig = sqrt(2 .* prop.L .* D ./ prop.v_bar); % diffusive spreading parameter

%-- Evaluate relevant functions ------------------------------------------%
rho_fun = @(G, r) (G - r) ./ (sqrt(2) .* sig); % recurring quantity
kap_fun = @(G, r) ...
    (G - r) .* erf(rho_fun(G, r)) + ...
    sig .* sqrt(2 / pi) .* exp(-rho_fun(G, r) .^ 2); % define function for kappa

%-- Evaluate the transfer function and its terms -------------------------%
%   real(...) is added for Case 2C.
K22 = kap_fun(real(G0(prop.r2)), prop.r2);
K21 = kap_fun(real(G0(prop.r2)), prop.r1);
K12 = kap_fun(real(G0(prop.r1)), prop.r2);
K11 = kap_fun(real(G0(prop.r1)), prop.r1);
Lambda = -1 ./ (4 .* prop.del) .* (K22 - K12 - K21 + K11);
Lambda(K22 > 1e2) = 0; % remove cases with large values out of error fun. eval.
Lambda(abs(Lambda) < 1e-10) = 0; % remove cases with roundoff error

end
