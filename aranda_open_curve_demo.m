%function aranda_open_curve_demo()
    % 1. Define an Open Curve (e.g., a noisy half-circle)
    t = linspace(0, pi, 50); 
    x = cos(t) + 0.02*randn(size(t)); % Add noise
    y = sin(t) + 0.02*randn(size(t)); 
    
    P_open = [x; y]; % 2xM matrix
    
    % 2. The "Staib & Duncan" Trick: Retrace to close the loop
    % We append the points in reverse order (excluding endpoints to avoid duplication)
    P_closed = [P_open, fliplr(P_open(:, 2:end-1))];
    
    % N is now the number of agents/points in the "virtual" closed formation
    N = size(P_closed, 2); 
    
    % 3. Calculate Aranda's Fourier Descriptors (Eq. 6)
    % D will be a 2xN matrix where column k is descriptor d_k
    D = zeros(2, N);
    
    for k = 0:N-1
        sum_vec = [0; 0];
        for i = 0:N-1
            % Aranda Eq (6): d_k = Sum( R(k*i) * p_i )
            % Note: Aranda indices are 0-based; Matlab is 1-based for storage
            p_i = P_closed(:, i+1);
            sum_vec = sum_vec + get_R(k*i, N) * p_i;
        end
        D(:, k+1) = sum_vec;
    end
    
    % 4. Filter: Keep only Low Frequencies (Aranda's Formation Specification)
    % Select H (number of harmonics). Aranda suggests odd numbers e.g., 3, 5.
    H = 3; 
    
    % Define the set of indices to KEEP (Aranda Eq. 10)
    % H_indices = {0, 1, ..., (H-1)/2,  N-(H-1)/2, ..., N-1}
    limit = (H-1)/2;
    keep_indices = [0:limit, (N-limit):N-1];
    
    % Create filtered descriptors (zero out high frequencies)
    D_filtered = zeros(2, N);
    % Map logic 0-based indices to 1-based Matlab indices
    D_filtered(:, keep_indices + 1) = D(:, keep_indices + 1);
    
    % 5. Reconstruct the Curve (Inverse DFT - Aranda Eq. 7)
    P_reconstructed = zeros(2, N);
    
    for i = 0:N-1
        sum_vec = [0; 0];
        for k = 0:N-1
            d_k = D_filtered(:, k+1);
            % IDFT uses R(k*i)' (Transpose of rotation matrix)
            sum_vec = sum_vec + get_R(k*i, N)' * d_k;
        end
        P_reconstructed(:, i+1) = (1/N) * sum_vec;
    end
    
    % 6. Extract the Open Curve
    % The first M points correspond to the original open segment
    M = size(P_open, 2);
    P_final_open = P_reconstructed(:, 1:M);

    % --- Visualization ---
    figure; hold on; axis equal; grid on;
    
    % Plot Original Noisy Points
    plot(P_open(1,:), P_open(2,:), 'k.', 'MarkerSize', 15, 'DisplayName', 'Original Points');
    
    % Plot the "Virtual" Closed Loop (for understanding)
    plot(P_reconstructed(1,:), P_reconstructed(2,:), 'g--', 'LineWidth', 1, 'DisplayName', 'Virtual Closed Loop');
    
    % Plot the Resulting Smooth Open Curve
    plot(P_final_open(1,:), P_final_open(2,:), 'r-', 'LineWidth', 2, 'DisplayName', 'Aranda Open Curve FD');
    
    legend;
    title(['Aranda FDs on Open Curve (H = ' num2str(H) ')']);
    xlabel('x'); ylabel('y');
%end

% Helper: Aranda's Rotation Matrix R(l) - Eq (5)
function R = get_R(l, N)
    % Aranda defines angles in radians: -2*pi*l / N
    theta = -2 * pi * l / N;
    R = [cos(theta), -sin(theta);
         sin(theta),  cos(theta)];
end