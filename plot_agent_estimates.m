function plot_agent_estimates(trajectory, time_points, params, agents_to_plot)
% PLOT_AGENT_ESTIMATES Plots True vs Estimated Fourier Descriptors
%
% Usage:
%   plot_agent_estimates(trajectory, time_points, params, [1, 3])
%
% Inputs:
%   trajectory      : [T x StateDim] matrix from simulate_rk4
%   time_points     : [T x 1] time vector
%   params          : struct containing N, H, MH, etc.
%   agents_to_plot  : vector of agent indices (1-based, e.g., [1, 2])

    N = params.N;
    H = params.H;
    MH = params.MH;
    xdim = params.xdim_agent;
    
    T_steps = length(time_points);
    num_agents_plot = length(agents_to_plot);

    % Pre-allocate storage for true descriptors
    d_true_hist = zeros(2*H, T_steps);
    
    % Pre-allocate storage for estimates: cell array of matrices
    % Each cell holds the history for one agent in 'agents_to_plot'
    d_est_hist = cell(1, num_agents_plot);
    for k = 1:num_agents_plot
        d_est_hist{k} = zeros(2*H, T_steps);
    end

    % --- Loop through trajectory to reconstruct data ---
    for t = 1:T_steps
        X = trajectory(t, :)';
        
        % 1. Extract Positions (p)
        p = zeros(2*N, 1);
        for i = 1:N
            base = (i-1)*xdim;
            p(2*i-1 : 2*i) = X(base+1 : base+2);
        end
        
        % 2. Calculate TRUE Global Descriptors (Ground Truth)
        % d = MH * p
        d_true = MH * p; 
        d_true_hist(:, t) = d_true;

        % 3. Extract Estimates for selected agents
        % Estimator state 'z' starts after all agent states
        z_start_idx = N * xdim; 
        z_all = X(z_start_idx+1 : end);
        Zmat = reshape(z_all, 2*H, N); % Column i is z_i
        
        for k = 1:num_agents_plot
            agent_idx = agents_to_plot(k);
            
            % Reconstruct local estimate: d_hat_i = z_i + r_i
            % r_i = N * MH(:, cols) * p_i
            
            p_i = p(2*agent_idx-1 : 2*agent_idx);
            cols = (2*agent_idx-1) : (2*agent_idx);
            
            r_i = N * (MH(:, cols) * p_i);
            z_i = Zmat(:, agent_idx);
            
            d_est_i = z_i + r_i;
            d_est_hist{k}(:, t) = d_est_i;
        end
    end

    % --- Plotting ---
    % We will create one figure per selected agent to avoid clutter
    
    % Define descriptor indices for convenience
    % Convention: DC=1:2, k=1 (d1)=3:4, k=-1 (d_1)=2H-1:2H (usually)
    % Based on 'build_low_freq_projection':
    % Indices correspond to frequencies: 0, 1, ..., (H-1)/2, -(H-1)/2, ... -1
    
    idx_d0 = 1:2;
    idx_d1 = 3:4;
    idx_dn1 = (2*H-1):(2*H); % Last 2 components are usually frequency -1

    for k = 1:num_agents_plot
        agent_idx = agents_to_plot(k);
        est_data = d_est_hist{k};
        
        figure('Color','w', 'Name', sprintf('Agent %d Estimation', agent_idx));
        
        % Subplot 1: Fundamental Ellipse (k=1) Norm
        subplot(3,1,1); hold on; grid on;
        norm_d1_true = vecnorm(d_true_hist(idx_d1, :));
        norm_d1_est  = vecnorm(est_data(idx_d1, :));
        plot(time_points, norm_d1_true, 'k-', 'LineWidth', 1.5);
        plot(time_points, norm_d1_est, 'r--', 'LineWidth', 1.5);
        ylabel('||d_1||');
        title(sprintf('Agent %d: Fundamental Component (k=1)', agent_idx));
        legend('True', 'Estimated');
        
        % Subplot 2: Centroid (k=0) Norm
        subplot(3,1,2); hold on; grid on;
        norm_d0_true = vecnorm(d_true_hist(idx_d0, :));
        norm_d0_est  = vecnorm(est_data(idx_d0, :));
        plot(time_points, norm_d0_true, 'k-', 'LineWidth', 1.5);
        plot(time_points, norm_d0_est, 'b--', 'LineWidth', 1.5);
        ylabel('||d_0||');
        title('Centroid Component (k=0)');
        
        % Subplot 3: Error Norm (Total Estimation Error)
        subplot(3,1,3); hold on; grid on;
        error_vec = d_true_hist - est_data;
        total_error = vecnorm(error_vec);
        plot(time_points, total_error, 'm-', 'LineWidth', 1.5);
        xlabel('Time (s)');
        ylabel('||d_{true} - d_{est}||');
        title('Total Estimation Error');
    end
end

