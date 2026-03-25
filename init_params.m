function params = init_params()

params.N = 7; % 7 for nice plot
params.H = 3;

params.xdim_agent = 3; % do not change this

params.vf = 0.5; % Agent linear velocity forward/backwards
params.vb = 0.5; 

params.omegal = pi/4; % Agent angular velocity left/right
params.omegar = pi/4;

params.kappa = 20 * params.N * params.vf; % estimator parameter

params.A_star = -100 ; % Desired area of the formation

end
