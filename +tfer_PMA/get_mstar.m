
% GET_MSTAR Get mass setpoint from provided speeds and voltages
% Author:   Timothy Sipkens, 2019-15-07
%=========================================================================%

function [m_star] = get_mstar(prop,V,omega1,omega2)

e = 1.60218e-19; % electron charge [C]

omega_hat = omega2./omega1;

alpha = omega1.*(prop.r_hat.^2-omega_hat)./(prop.r_hat.^2-1);
beta = omega1.*prop.r1.^2.*(omega_hat-1)./(prop.r_hat^2-1);

m_star = V./...
    (log(1/prop.r_hat)./e.*(alpha.*prop.rc+beta./prop.rc).^2);

end

