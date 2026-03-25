function plot_fd_area_metrics(traj, tspan, params)
% plot_fd_area_metrics
% Plots ||d_{-1}||, ||d_0||, ||d_1|| and the formation area over time
%
% Inputs:
%   traj   : [T x nX] state trajectory
%   tspan  : [T x 1] time vector
%   params : struct with fields
%       N, H, xdim_agent, MH
%       (optional) A_star

    N    = params.N;
    H    = params.H;
    xdim = params.xdim_agent;
    MH   = params.MH;

    assert(mod(H,2)==1, 'H must be odd.');

    T = length(tspan);

    % Storage
    d0_norm  = zeros(T,1);
    d1_norm  = zeros(T,1);
    d_1_norm = zeros(T,1);
    area_val = zeros(T,1);

    K = (H-1)/2;

    for t = 1:T
        X = traj(t,:)';

        %% 1. Extract agent positions p
        p = zeros(2*N,1);
        for i = 1:N
            base = (i-1)*xdim;
            p(2*i-1:2*i) = X(base+1:base+2);
        end

        %% 2. Extract estimator states z
        z_all = X(N*xdim+1:end);
        Zmat  = reshape(z_all, 2*H, N);

        %% 3. Reconstruct local descriptor estimates
        R = zeros(2*H, N);
        for i = 1:N
            cols = (2*i-1):(2*i);
            R(:,i) = N * (MH(:,cols) * p(2*i-1:2*i));
        end
        Dhat = Zmat + R;

        %% 4. Consensus descriptor (global estimate)
        d = mean(Dhat, 2);

        %% 5. Extract low-frequency descriptors
        % Indexing convention:
        % [ d0 | d1 | d2 ... dK | d_{N-K} ... d_{N-1} ]

        d0   = d(1:2);
        d1   = d(3:4);
        d_1  = d(end-1:end);

        d0_norm(t)  = norm(d0);
        d1_norm(t)  = norm(d1);
        d_1_norm(t) = norm(d_1);

        %% 6. Area computation (general odd-H case)
        A_sum = 0;
        for k = 1:K
            idx_pos = 2*k + (1:2);
            idx_neg = 2*H - 2*k + (1:2);

            dk_pos = d(idx_pos);
            dk_neg = d(idx_neg);

            A_sum = A_sum + k * (norm(dk_pos)^2 - norm(dk_neg)^2);
        end

        area_val(t) = (pi / N^2) * A_sum;
    end

    %% === Plotting ===
    figure('Color','w','Name','Fourier Descriptor Metrics');

    % --- Descriptor norms ---
    subplot(2,1,1); hold on;
    plot(tspan, d_1_norm, 'LineWidth', 1.5);
    plot(tspan, d0_norm,  'LineWidth', 1.5);
    plot(tspan, d1_norm,  'LineWidth', 1.5);
    grid on;
    xlabel('Time (s)');
    ylabel('Descriptor norm');
    title('Norm of k=1 (H=3 ellipse) Fourier descriptors');
    legend({'||d_{-1}||','||d_0||','||d_1||'}, 'Location','best');
    hold off;

    % --- Area ---
    subplot(2,1,2); hold on;
    plot(tspan, area_val, 'LineWidth', 2);
    if isfield(params,'A_star')
        yline(params.A_star,'--r','Target Area (A^*)','LineWidth',1.5, 'LabelHorizontalAlignment', 'left');
        legend({'A(t)','A^*'}, 'Location','best');
    else
        legend({'A(t)'}, 'Location','best');
    end
    
    % Plot reference area if available
    %if isfield(params, 'A_star')
    %    yline(params.A_star, '--r', 'Target Area (A^*)', 'LineWidth', 1.5, 'LabelHorizontalAlignment', 'left');
    %end

    grid on;
    xlabel('Time (s)');
    ylabel('Area');
    title('Formation area');
    hold off;
end
