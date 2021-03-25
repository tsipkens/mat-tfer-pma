
% MASSMOB  Fill in mass-mobility information.
%  Includes computing prop.m0, which is used for mp2zp.
%  
%  AUTHOR: Timothy Sipkens, 2021-03-25

function prop = prop_massmob(prop)

if ~isfield(prop, 'm0')
    if isfield(prop, 'm100')
        prop.m0 = prop.m100 * (1 / 100) ^ prop.Dm;
        
    elseif isfield(prop, 'rho0')
        prop.m0 = prop.rho0 * pi / 6 * 1e-27;
        
    elseif isfield(prop, 'rho100')
        prop.m100 = prop.rho100 * pi / 6 * (100e-9) ^ 3;
        prop.m0 = prop.m100 * (1 / 100) ^ prop.Dm;
        
    else
        error('Could not compute prop.m0.');
    end
end

prop.rho0 = prop.m0 * 6 / pi * 1e27;
prop.m100 = prop.m0 / (1 / 100) ^ prop.Dm;
prop.rho100 = prop.m100 * 6 / pi / (100e-9 ^ 3);

end
