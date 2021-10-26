
% TFER_2C_DIFF  Evaluates the transfer function for a PMA in Case D (w/ diffusion).
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

function [Lambda,G0] = tfer_2C_diff(sp,m,d,z,prop)

[~,~,D] = parse_inputs(sp,m,d,z,prop); % get diffusion coeff.
sig = sqrt(2.*prop.L.*D./prop.v_bar); % diffusive spreading parameter


%-- Evaluate relevant functions ------------------------------------------%
[~,G0] = tfer_2C(sp,m,d,z,prop);
    % get G0 function for this case

Lambda = tfer_diff(G0, D, prop);  % apply general diffusive form

end
