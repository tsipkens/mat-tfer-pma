
% GET_MSTAR Get mass setpoint for a CPMA from provided speeds and voltages.
% Author:   Timothy Sipkens, 2019-15-07
%=========================================================================%

function [m_star,prop] = get_mstar(prop,V,omega1,omega2)
%-------------------------------------------------------------------------%
% Inputs:
%   prop    Struct containing physical dimensions of CPMA
%   V       Operating voltage of the CPMA [V]
%   omega1  Rotational speed of inner electrode [rad/s]
%   omega2  Rotation speed of outer electrode [rad/s]
%
% Ouputs:
%   m_star  Mass setpoint [kg]
%   prop    Updated struct containing physical dimensions of CPMA, where
%           omega_hat is updated
%-------------------------------------------------------------------------%

e = 1.60218e-19; % electron charge [C]

omega_hat = omega2./omega1;
prop.omega_hat = omega_hat; % update prop to reflect setpoints

alpha = omega1.*(prop.r_hat.^2-omega_hat)./(prop.r_hat.^2-1);
beta = omega1.*prop.r1.^2.*(omega_hat-1)./(prop.r_hat^2-1);

m_star = V./...
    (log(1/prop.r_hat)./e.*(alpha.*prop.rc+beta./prop.rc).^2);

end

