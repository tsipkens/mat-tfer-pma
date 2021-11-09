
% PROP_UPDATE_MASSMOB  A simple utility to update the mass-mobility parameters.
%  
%  PROP = prop_update_massmob(PROP, F1, V1, F2, V2) takes an input
%  property structure and replace the mass-mobility relation parameters
%  with the susbequent inputs. F* variables are strings containing the
%  field names, while V* contains the corresponds values for those 
%  parameters. For example, F1 = 'Dm' updates the mass-mobility exponent to
%  the value given as V1. A pair of mass-mobility parameters is required to
%  constrain the relation. 
%  
%  NOTE: The function simply removes the existing parameters and call
%  prop_massmod(PROP) to reintroduce the parameters. See that function for
%  further details. 
%  
%  ------------------------------------------------------------------------
%  
%  AUTHOR: Timothy Sipkens, 2021-10-06

function prop = prop_update_massmob(prop, f1, v1, f2, v2)

% Remove the relevant fields, to be replaced.
prop = rmfield(prop, {'Dm', 'm0', 'rho0', 'm100', 'rho100'});

% Add new values.
prop.(f1) = v1;
prop.(f2) = v2;

% Fill out remaining values.
prop = prop_massmob(prop);

end

