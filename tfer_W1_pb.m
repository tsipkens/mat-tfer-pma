
% TFER_W1_PB    Evaluates the transfer function for a PMA in Case E (w/ parabolic flow).
% Author:       Timothy Sipkens, 2019-03-21
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

function [Lambda, G0] = tfer_W1_pb(sp, m, d, z, prop)

[tau, C0, ~, rs] = parse_inputs(sp, m, d, z, prop);
        % parse inputs for common parameters

%-- Estimate recurring quantities ----------------------------------------%
A1 = -3 .* prop.v_bar ./ (4 .* tau .* [sp.omega1]' .^ 4 .* prop.del .^ 2);
A2 = [sp.omega1]' .^ 2 .* (prop.rc .^ 2 - prop.del .^ 2) + C0 ./ m;
A3 = 4 .* prop.rc .* [sp.omega1]' .* sqrt(C0 ./ m);


% Loop over setpoints.
Lambda = zeros(size(A1));
for jj=1:size(Lambda,1)
    
    A4 = @(r,ii) sp(jj).omega1 .^ 2 .* (r .^ 2 - 4 .* prop.rc .* r);
    
    %-- Set up F function for minimization -------------------------------%
    F = @(r,ii) A1(jj,ii) .* (A2(jj,ii) .* log(C0(jj) ./ m(ii) - ...
        (sp(jj).omega1 .* r) .^2 )+...
        A3(jj,ii) .* atanh(sqrt(m(ii) ./ C0(jj)) .* sp(jj).omega1 .* r) + ...
        A4(r, ii));
    min_fun = @(rL,r0,ii) F(rL,ii) - F(r0,ii) - prop.L;

    %-- Evaluate G0 and transfer function --------------------------------%
    G0 = @(r) G_fun(min_fun, r, rs(jj,:), ...
        prop.r1, prop.r2, sp(jj).alpha, sp(jj).beta);

    ra = min(prop.r2, max(prop.r1, G0(prop.r1)));
    rb = min(prop.r2, max(prop.r1, G0(prop.r2)));

    Lambda(jj,:) = 3/4 .* (rb - ra) ./ prop.del - ...
        1/4 .* ((rb - prop.rc) ./ prop.del) .^ 3 + ...
        1/4 .* ((ra - prop.rc) ./ prop.del) .^ 3;
end

end

