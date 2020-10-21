
% MP2ZP     Calculate electric mobility from a vector of particle mass.
% Author:   Timothy Sipkens, 2019-01-02
% 
% Inputs:
%   m           Particle mass
%   z           Integer charge state
%   T           System temperature
%   P           System pressure
%   prop        CPMA/DMA properties structure
%
% Outputs:
%   B           Mechanical mobility
%   Zp          Electromobility
%   d           Mobility diameter (implied by mass-mobility relation)
% 
% Note:
%   Uses mass-mobility relationship to first convert to a mobility
%   diameter and then estimates the mobility using dm2zp.
%=========================================================================%

function [B,Zp,d] = mp2zp(m,z,T,P,prop)


%-- Parse inputs ---------------------------------------------------------%
if ~exist('T', 'var'); T = []; end
if ~exist('P', 'var'); P = []; end

if ~exist('prop', 'var'); prop = []; end
if or(isempty(prop),...
        ~and(isfield(prop,'m0'),...
        isfield(prop,'Dm'))) % get parameters for the mass-mobility relation
    error(['Please specify the mass-mobility relation parameters ',...
        'in the prop structure.']);
end
%-------------------------------------------------------------------------%


d = 1e-9 .* (m ./ prop.m0) .^ (1 / prop.Dm);
    % use mass-mobility relationship to get mobility diameter

    
%-- Use mobility diameter to get particle electro and mechanical mobl. ---%
if or(isempty(T),isempty(P))
    [B, Zp] = dm2zp(d, z);
else
    [B, Zp] = dm2zp(d, z, T, P);
end

end

