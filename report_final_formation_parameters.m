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
