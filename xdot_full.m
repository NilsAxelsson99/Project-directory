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

% compute local reference r_i for each agent (2H x 1)
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
    
% comment our this to get old file
% Applies the minimal-change projection to reach params.A_star
%         if isfield(params, 'A_star') && H >= 3
%             % Extract descriptors: d1 is at index 3:4, dN-1 is at the end
%             d1  = dhat_i(3:4);       
%             dN1 = dhat_i(end-1:end); 
% 
%             % Calculate current local area estimate: A = (pi/N^2)(||d1||^2 - ||dN-1||^2)
%             A_hat = (pi / N^2) * (norm(d1)^2 - norm(dN1)^2);
% 
%             % Compute correction factor alpha based on the linearized gradient
%             % denominator is proportional to the 'energy' of the shape descriptors
%             denom = (2 * pi / N^2) * (norm(d1)^2 + norm(dN1)^2);
%             alpha = (params.A_star - A_hat) / (denom + 1e-6); % Add eps for stability
% 
%             % Update descriptors locally using radial scaling
%             dhat_i(3:4)       = d1  * (1 + alpha);
%             dhat_i(end-1:end) = dN1 * (1 - alpha);
%         end

        % --------------------------------------

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
       % alpha = 1 / (denom); 

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
