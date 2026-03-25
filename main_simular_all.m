%main_simulation.m  
% │          
% │
% ├── dynamics/
% │   ├── xdot_full.m
% │   ├── agent_unicycle.m (NOT HERE)
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
params.centroid_ref = centroid_reference("lemniscate");

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
plot_fd_area_metrics(trajectory, time_points, params)
plot_trajectories(trajectory, params);
animate_simulation(trajectory, time_points, params);
report_final_formation_parameters(trajectory, params)

plot_dH_agent(trajectory, time_points, params, 3) % test

%% Main parts to make everything run
function x_next = rk4_step(f, t, x, h, params)
    k1 = f(t, x, params);
    k2 = f(t + h/2, x + (h/2) * k1, params);
    k3 = f(t + h/2, x + (h/2) * k2, params);
    k4 = f(t + h,   x + h * k3, params);
    x_next = x + (h/6) * (k1 + 2*k2 + 2*k3 + k4);
end

function [trajectory, time_points] = simulate_rk4(f, X0, t0, T, h, params)

% ----- Setup -------- %
num_steps = round((T - t0)/h);
state_dim = length(X0);

trajectory = zeros(num_steps+1, state_dim);
time_points = zeros(num_steps+1,1);

trajectory(1,:) = X0';
time_points(1) = t0;

X = X0; t = t0;

% ----- Integrate over all time ----- %
    for k = 1:num_steps 
        X = rk4_step(f, t, X, h, params);
        t = t + h;
        trajectory(k+1,:) = X';
        time_points(k+1) = t;
    end

end

function dX = xdot_full(t, X, params)
% dX = xdot_full(t, X, params)
% X = [agents; z_all]
N = params.N;
H = params.H;
xdim_agent = params.xdim_agent;

MH = params.MH;       % 2H x 2N
Neighbors = params.Neighbors;
kappa = params.kappa;

% saturation velocities
vf = params.vf; 
vb = params.vb;
omegal = params.omegal;
omegar = params.omegar;

% unpack agent states into p (2N) and thetas (N)
p = zeros(2*N,1);
thetas = zeros(N,1);
for i = 1:N
    base = (i-1)*xdim_agent;
    p(2*i-1) = X(base + 1);
    p(2*i)   = X(base + 2);
    thetas(i) = X(base + 3);
end

% unpack z_all
z_all = X(N*xdim_agent + 1 : end);   % length N*2H
Zmat = reshape(z_all, 2*H, N);

% compute local r_i for each agent (2H x 1)
R = zeros(2*H, N);
for i = 1:N
    p_i = p(2*i-1:2*i);
    cols = (2*i-1):(2*i);
    R(:, i) = N * (MH(:, cols) * p_i);
end

% local raw estimates dhat_i = z_i + r_i
Dhat = Zmat + R;   % 2H x N

% estimator dynamics (can use smoothed sign/tanh)
Zdot = zeros(2*H, N);
for i = 1:N
    sumterm = zeros(2*H,1);
    for j = Neighbors{i}
        delta = Dhat(:, j) - Dhat(:, i);
        sumterm = sumterm + sign( delta ); % use smooth approx of sign?
    end
    Zdot(:, i) = kappa * sumterm;
end

% Now compute agent control using each agent's local d_H^{(i)},
% but override the DC components (1:2) with centroid_ref(t)
agent_dot = zeros(N * xdim_agent, 1);

% Loop over all agents to calculate their respective xdot
for i = 1:N 
    dhat_i = Dhat(:, i);

%   ---  Generalized area-regulation for arbitrary odd H ---
    if isfield(params, 'A_star') && H >= 3
        K = (H - 1) / 2;
        A_hat_sum = 0;
        grad_norm_sum = 0;
    
        % 1) Calculate current area and the gradient normalization factor
        for k = 1:K
            % Mapping indices based on build_low_freq_projection order:
            % [DC, k=1, k=2... k=K, N-K, N-K+1... N-1]
            idx_pos = 2*k + (1:2);
            idx_neg = 2*H - 2*k + (1:2); % N-k is at the end of the vector
    
            dk_pos = dhat_i(idx_pos);
            dk_neg = dhat_i(idx_neg);
    
            % Area contribution: k*(||dk_pos||^2 - ||dk_neg||^2)
            A_hat_sum = A_hat_sum + k * ( norm(dk_pos)^2 - norm(dk_neg)^2 );
    
            % Gradient energy/2-norm: k^2*(||dk_pos||^2 + ||dk_neg||^2)
            % This forms the denominator for the projection step
            grad_norm_sum = grad_norm_sum + k^2 * (norm(dk_pos)^2 + norm(dk_neg)^2);
        end
    
        % Area computation
        A_hat = (pi / N^2) * A_hat_sum;
    
        % Denominator based on the norm of the total gradient vector
        denom = (2 * pi / N^2) * grad_norm_sum; % no square since it cancels the numerators similar term
        
        % Local area error for each agent
        eA = params.A_star - A_hat;
    
        % Compute correction factor alpha = (A_d - A(d) / (norm{nabla A(d}^2)
        %alpha = (params.A_star - A_hat) / (denom + 1e-6); 
        alpha = 1 / (denom + 1e-6); 
    
        % 2) Update descriptors locally using the frequency-weighted gradient
        for k = 1:K
            idx_pos = 2*k + (1:2);
            idx_neg = 2*H - 2*k + (1:2);
    
            % Move in the direction of the gradient to increase/decrease area
            dhat_i(idx_pos) = dhat_i(idx_pos) * (1 + k * eA * alpha);
            dhat_i(idx_neg) = dhat_i(idx_neg) * (1 - k * eA * alpha);
        end
    end
    %---------------------------------------------------------------------
    

    % override DC with commanded centroid to include tracking
    c_ref = params.centroid_ref(t);
    dhat_i(1:2) = N * c_ref;   % Correct scaling for DC mode, otherwise all agents go to, e.g., 10/N when ref is 10 if N agents
   

    % reconstruct projected positions from this (1/N * MH' * dhat_i)
    pH_hat = (1/N) * (MH' * dhat_i);   % 2N x 1
    pH_i = pH_hat(2*i-1:2*i);
    
    p_i = p(2*i-1:2*i);
    f_i = pH_i - p_i;   % single-integrator desired velocity (2x1) (=error)

    th = thetas(i);
    e_long = [cos(th), sin(th)] * f_i;     % scalar
    e_ang  = [-sin(th), cos(th)] * f_i;    % scalar

    v     = sat_v(e_long, vf, vb);
    omega = sat_omega(e_ang, omegal, omegar);

    base = (i-1)*xdim_agent;
    agent_dot(base + 1) = v * cos(th);
    agent_dot(base + 2) = v * sin(th);
    agent_dot(base + 3) = omega;
end

% pack derivatives: [agent_dot; zdot_flat]
dX = [agent_dot; reshape(Zdot, [], 1)];
end

% --- saturation functions ---
function v_s = sat_v(v, vf, vb)
    v_s = min(max(v, -vb), vf);
end

function w_s = sat_omega(w, omegal, omegar)
    w_s = min(max(w, -omegar), omegal);
end

%% SUPPORT FUNCTIONS
function report_final_formation_parameters(traj, params)
% report_final_formation_parameters
% Extracts the final consensus Fourier descriptors, computes the
% formation area (general odd-H case), and reports the k=1 ellipse
% parameters following Sjöstrand's parametrization.

    %% Parameters
    N    = params.N;
    H    = params.H;
    xdim = params.xdim_agent;
    MH   = params.MH;

    assert(mod(H,2)==1, 'H must be odd.');
    K = (H-1)/2;

    %% Final state
    Xf = traj(end,:)';

    %% 1. Extract agent positions
    p = zeros(2*N,1);
    for i = 1:N
        base = (i-1)*xdim;
        p(2*i-1:2*i) = Xf(base+1:base+2);
    end

    %% 2. Extract estimator states
    z_all = Xf(N*xdim+1:end);
    Zmat  = reshape(z_all, 2*H, N);

    %% 3. Reconstruct local descriptor estimates
    R = zeros(2*H, N);
    for i = 1:N
        cols = (2*i-1):(2*i);
        R(:,i) = N * (MH(:,cols) * p(2*i-1:2*i));
    end

    Dhat = Zmat + R;

    %% 4. Consensus descriptor vector
    d = mean(Dhat, 2);

    %% 5. DC component (centroid)
    d0 = d(1:2);
    c  = d0 / N;

    %% 6. Area computation (general odd-H)
    Area = 0;
    for k = 1:K
        idx_pos = 2*k + (1:2);
        idx_neg = 2*H - 2*k + (1:2);

        dk_pos = d(idx_pos);
        dk_neg = d(idx_neg);

        Area = Area + (pi/N^2) * k * (norm(dk_pos)^2 - norm(dk_neg)^2);
    end

    %% 7. k = 1 descriptors (ellipse parameters)
    d1  = d(3:4);
    d_1 = d(end-1:end);

    % Sjöstrand / Aranda mapping
    a1 = (d1(1) + d_1(1)) / N;
    b1 = (d_1(2) - d1(2)) / N;
    c1 = (d1(2) + d_1(2)) / N;
    d1c = (d1(1) - d_1(1)) / N;

    %% 8. Reporting
    sep = [newline repmat('=',1,50) newline];
    fprintf(sep);
    fprintf('FINAL FORMATION SUMMARY (H = %d)\n', H);
    fprintf(sep);

    fprintf('Consensus Fourier Descriptors:\n');
    fprintf('  d0  = [%.4f, %.4f]\n', d0(1), d0(2));
    fprintf('  d1  = [%.4f, %.4f]\n', d1(1), d1(2));
    fprintf('  d-1 = [%.4f, %.4f]\n', d_1(1), d_1(2));

    fprintf(' ||d1||  = %.4f\n', norm(d1));
    fprintf(' ||d-1|| = %.4f\n', norm(d_1));
    fprintf(' ||d0||  = %.4f\n', norm(d0));
    fprintf(' (H=3) Area = pi/(N^2) *( ||d1||^2-||d-1||^2)  = %.4f \n', ( pi/(N^2) ) * ( (norm(d1))^2 - (norm(d_1))^2 )  );

    fprintf('\nCentroid (d0/N) ):\n');
    fprintf('  c = [%.4f, %.4f]\n', c(1), c(2));

    fprintf('\nFundamental Ellipse (k = 1):\n');
    fprintf('  P = [ a1  b1 ; c1  d1 ]\n');
    fprintf('      [ %8.4f  %8.4f ]\n', a1, b1);
    fprintf('      [ %8.4f  %8.4f ]\n', c1, d1c);
    fprintf(' Area_check: A = %8.4f \n', pi*( (a1*d1c)-(c1*b1) ) );

    fprintf('\nArea:\n');
    if isfield(params,'A_star')
        fprintf('  Target Area A*     = %.4f\n', params.A_star);
    end
    fprintf('  Computed Area A    = %.4f\n', Area);

    fprintf(sep);
end

function plot_dH_agent(traj, tspan, params, agent_id)
% plot_dH_agent
% ------------------------------------------------------------
% Plots the local descriptor estimate dH^(i)(t) for a selected agent
%
% INPUTS:
%   traj      : [T x nX] state trajectory
%   tspan     : [T x 1] time vector
%   params    : struct with fields N, H, xdim_agent, MH
%   agent_id  : index of agent to visualize
% ------------------------------------------------------------

N    = params.N;
H    = params.H;
xdim = params.xdim_agent;
MH   = params.MH;

assert(agent_id >= 1 && agent_id <= N, ...
    'agent_id must be between 1 and %d', N);

Tsteps = length(tspan);

% Storage for selected agent descriptor
dH_agent = zeros(2*H, Tsteps);

for t = 1:Tsteps

    X = traj(t,:)';

    % ---- Extract positions ----
    p = zeros(2*N,1);
    for i = 1:N
        base = (i-1)*xdim;
        p(2*i-1:2*i) = X(base+1:base+2);
    end

    % ---- Extract estimator states ----
    z_all = X(N*xdim+1:end);
    Zmat  = reshape(z_all, 2*H, N);

    % ---- Compute local reference R ----
    R = zeros(2*H, N);
    for i = 1:N
        cols = (2*i-1):(2*i);
        R(:,i) = N * (MH(:,cols) * p(2*i-1:2*i));
    end

    % ---- Local descriptor estimates ----
    Dhat = Zmat + R;

    % Store chosen agent
    dH_agent(:,t) = Dhat(:,agent_id);
end

%% ===== Plotting =====
figure('Color','w');
plot(tspan, dH_agent.', 'LineWidth', 1.2);
grid on;

xlabel('Time (s)');
ylabel('Descriptor components');
title(sprintf('Local Descriptor Estimate d_H^{(%d)}(t)', agent_id));

legend(arrayfun(@(k) sprintf('d_{%d}',k), ...
        1:2*H, 'UniformOutput', false), ...
        'Location','bestoutside');

end


