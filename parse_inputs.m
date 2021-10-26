
% PARSE_INPUTS A function to evaluate setpoint parameters including C0, alpha, and beta.
% Author:  Timothy A. Sipkens, 2019-05-01
% 
% Required variables:
%   sp      Mass corresponding to the measurement set point of the APM
%   d           Struct containing mutliple setpoint parameters (V, alpha, etc.)
%   m           Particle mass, can be vector of same length as d
%   z           Integer charge state, scalar
%   prop        Properties of particle mass analyzer
%
% Sample outputs:
%   C0          Summary parameter for the electrostatic force
%   tau         Product of mechanical mobility and particle mass
%   D           Diffusion coefficient for specified particles
%   rs          Equilibrium radius
% 
% Notes:    As a script, this code uses variables currently in the 
%           workspace. This script is also used to parse some of the inputs 
%           to the various transfer functions, including the existence of 
%           the integer charge state and particle mobility. 
%=========================================================================%

function [tau, C0, D, rs] = parse_inputs(sp, m, d, z, prop)


%-- Parse inputs ---------------------------------------------------------%
if isempty(z); z = 1; end % if integer charge is not specified, use z = 1


%-- Set up mobility calculations -----------------------------------------%
e = 1.60218e-19; % electron charge [C]
q = z .* e; % particle charge


%-- Evaluate mechanical mobility -----------------------------------------%
if isempty(d) % if mobility diameter is NOT specified
    warning('Invoking mass-mobility relation to determine Zp.');
    B = mp2zp(m, z, prop.T, prop.p, prop);
else % if mobility diameter is specified
    B = dm2zp(d, z, prop.T, prop.p);
end


%-- Evaluate output parameters -------------------------------------------%
tau = B .* m;
D = prop.D(B) .* z; % diffusion as a function of mechanical mobiltiy and charge state
C0 = [sp.V]' .* q ./ log(1/prop.r_hat); % calcualte recurring C0 parameter

if nargout>=4 % if required, calculate equilbirium radius
% Note: Whether to pick the +ive of -ive root for rs is chosen based on a
% heuristic approach. Specifically, the root closer to the centerline is
% chosen, except when the -ive root is zero (which is the case for APM 
% conditions, where the +ive root should always be used).
    
    % evaluate +ive and -ive roots
    r_m = (sqrt(C0 ./ m) - ...
        sqrt(C0 ./ m - 4 .* [sp.alpha]' .* [sp.beta]')) ./ (2 .* [sp.alpha]');
    r_p = (sqrt(C0 ./ m) + ...
        sqrt(C0 ./ m - 4 .* [sp.alpha]' .* [sp.beta]')) ./ (2 .* [sp.alpha]');
    
    % determine which root is closer to centerline radius
    bo = abs(r_m-prop.rc) > abs(r_p-prop.rc);
    bo(r_m==0) = 1; % avoid zero values for APM case
    
    % assign one of the two roots to rs
    rs = r_m; % by default use -ive root
    rs(bo) = r_p(bo); % if closer to +ive root, use +ive root
    
    
    % zero out cases where no equilibrium radius exists (also removes complex numbers)
    rs(C0 ./ m < (4 .* [sp.alpha]' .* [sp.beta]')) = 0;
end


end
