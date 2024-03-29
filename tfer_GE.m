
% TFER_GE   Evaluates the transfer function for a PMA in Case F.
% Author:   Timothy Sipkens, 2019-03-21
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

function [Lambda,G0] = tfer_GE(sp, m, d, z, prop)

[tau, C0, ~, rs] = parse_inputs(sp, m, d, z, prop);
        % parse inputs for common parameters

%-- Estimate recurring quantities ----------------------------------------%
C6 = 2 .* [sp.alpha]' .* [sp.beta]' - C0 ./ m;
C7 = sqrt(4 .* [sp.alpha]' .^ 2 .* [sp.beta]' .^ 2 - C6 .^ 2);

A1 = prop.v_bar ./ (4 .* tau .* [sp.alpha]' .^ 2);
A2 = 2 .* C6 ./ C7;


Lambda = zeros(size(A1));
for jj=1:size(Lambda,1)
    %-- Set up F function for minimization -------------------------------%
    F = @(r,ii) A1(jj,ii) .* (log(sp(jj).alpha .^ 2 .* r .^ 4 + ...
        sp(jj).beta .^ 2 + C6(jj,ii) .* r .^ 2) - ...
        A2(jj,ii) .* atan((2 .* sp(jj).alpha .^ 2 .* r .^ 2 + ...
        C6(jj,ii)) ./ C7(jj,ii)));
    min_fun = @(rL,r0,ii) F(rL,ii) - F(r0,ii) - prop.L;


    %-- Evaluate G0 and transfer function --------------------------------%
    G0 = @(r) G_fun(min_fun, r, rs(jj,:), ...
        prop.r1, prop.r2, sp(jj).alpha, sp(jj).beta);
    
    ra = min(prop.r2, max(prop.r1, G0(prop.r1)));
    rb = min(prop.r2, max(prop.r1, G0(prop.r2)));

    Lambda(jj,:) = (1 / (2 .* prop.del)).*(rb - ra);
end

end

