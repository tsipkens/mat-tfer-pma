
% PROP_UPDATE_FLOW  A simple utility to update the flow rate and v_bar.
%  Acts as a reminder that these two parameters are linked.
%  
%  AUTHOR: Timothy Sipkens, 2021-03-26

function prop = prop_update_flow(prop, Q)

prop.Q = Q;
prop.v_bar = prop.Q / prop.A;

end

