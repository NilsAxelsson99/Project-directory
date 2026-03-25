%main_simulation.m  
% │          
% │
% ├── dynamics/
% │   ├── xdot_full.m
% │   ├── agent_unicycle.m
% │
% ├── integration/
% │   ├───rk4_step.m
% │   └── simulate_rk4.m
% │
% ├── descriptors/
% │   ├── calc_M_real.m
% │   └── build_low_freq_projection.m
% │
% ├── graph/
% │   └── build_ring_graph.m
% │
% ├── initialization/
% │   ├── init_agents.m
% │   ├── init_estimators.m
% │   └── init_params.m
% │
% ├── plotting/
% │   ├── plot_descriptors.m
% │   ├── plot_trajectories.m
% │   └── animate_simulation.m
% │
% └── utils/
%     └── stack_positions.m
close all; clear; 
%clc;
rng(1); % Set seed for randomized agent starting positions

agents_to_plot = [2,6]; % agents to look at

%% Time
t0 = 0.0;  % start
T = 150.0; % stop
h = 0.01;  % timestep

%% Parameters
params = init_params();
params.h = h;

% -------- Reference trajectory to track: ---------------%
% Trajectory options: 
    % "ellipse"
    % "circle"
    % "static" % ("static",[5;5])
    % "line"
    % "lemniscate" = \infty sign
%params.centroid_ref = centroid_reference("static",[0;0]);
params.centroid_ref = centroid_reference("ellipse");

% Initial conditions
[x0_agent, p0] = init_agents(params);
z0 = init_estimators(params);
X0 = [x0_agent; z0];

% Descriptor matrices, only need MH at the moment
 [params.MH, params.invMH, params.M, params.invM, params.low_idx_0] = ...
     build_low_freq_projection(params.N, params.H);

% Communication graph
params.Neighbors = build_ring_graph(params.N);

% Integrate
[trajectory, time_points] = simulate_rk4( @xdot_full, X0, t0, T, h, params);

disp("Simulation complete");

%% Plotting
plot_descriptors(p0, params);
plot_fd_area_metrics(trajectory, time_points, params)
plot_trajectories(trajectory, params);
animate_simulation(trajectory, time_points, params);
report_final_formation_parameters(trajectory, params)
%plot_agent_estimates(trajectory, time_points, params, agents_to_plot)

