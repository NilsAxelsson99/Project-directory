% NOT IN USE AT THE MOMENT
% function xdot = agent_unicycle(i, x_agents, z_all, t, params)
% 
% base  = (i-1)*params.xdim_agent;
% x     = x_agents(base+1);
% y     = x_agents(base+2);
% theta = x_agents(base+3);
% 
% % reference centroid
% c_ref = params.centroid_ref(t);
% 
% % simple centroid tracking law (unchanged logic)
% u     = c_ref - [x; y];
% v     = cos(theta)*u(1) + sin(theta)*u(2);
% omega = -sin(theta)*u(1) + cos(theta)*u(2);
% 
% % saturation
% [v, omega] = saturation(v, omega, params);
% 
% % Agent dynamics (for integration by RK4)
% xdot = [v*cos(theta);
%         v*sin(theta);
%         omega];
% end
