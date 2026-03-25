function [W_A, w] = calculate_WA_HA(N, H)
    % calculate_area_weights - Computes the signed frequency weights and 
    % the Hessian matrix for the area gradient control term.
    %
    % Inputs:
    %   N - Total number of agents
    %   H - Number of low-frequency descriptors (must be odd)
    %
    % Outputs:
    %   W_A - [2H x 2H] block-diagonal Hessian matrix
    %   w   - [H x 1] vector of signed weights w_k
    
    if mod(H, 2) == 0
        error('H must be odd to include symmetric frequencies and the DC component.');
    end
    
    % Number of positive/negative frequency pairs
    num_freqs = (H - 1) / 2;
    
    % Define the set of frequency indices k in H
    % Standard order: [0, 1, ..., num_freqs, N-num_freqs, ..., N-1]
    k_indices = [0, 1:num_freqs, (N-num_freqs):(N-1)];
    
    % Initialize weights and matrix
    w = zeros(H, 1);
    W_A = zeros(2*H, 2*H);
    
    for i = 1:H
        k = k_indices(i);
        
        % Calculate signed weight w_k
        if k == 0
            w_k = 0;
        elseif k <= num_freqs
            w_k = k;         % Positive frequencies
        else
            w_k = k - N;     % Negative frequencies (e.g., N-1 becomes -1)
        end
        
        w(i) = w_k;
        
        % Fill the 2x2 block in the Hessian matrix
        % Gradient factor: (2 * pi * w_k) / N^2
        block_val = (2 * pi * w_k) / (N^2);
        
        row_idx = 2*i - 1;
        W_A(row_idx : row_idx+1, row_idx : row_idx+1) = block_val * eye(2);
    end
end