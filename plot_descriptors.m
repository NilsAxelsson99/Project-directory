function plot_descriptors(p0, params)
% Plot original constellation, reconstructed constellation,
% and Fourier descriptor magnitudes

N = params.N;
H = params.H;

M     = params.M;
MH    = params.MH;
invM  = params.invM;
low_idx_0 = params.low_idx_0;

freq_axis   = 0:N-1;
freq_axis_H = low_idx_0;

% Stack initial positions
p0_stack = reshape(p0, [2*N, 1]);

% Full descriptors
d0  = M * p0_stack;
d0x = d0(1:2:end);
d0y = d0(2:2:end);

% Low-frequency descriptors
d0H  = MH * p0_stack;
d0Hx = d0H(1:2:end);
d0Hy = d0H(2:2:end);

% Reconstructed positions
p0_rec  = invM * d0;
p0_recx = p0_rec(1:2:end);
p0_recy = p0_rec(2:2:end);

%% Figure 1
figure(1);

    subplot(1,3,1)
    plot(p0(1,:), p0(2,:), 'ko', 'LineWidth',1.5); hold on;
    plot(p0_recx, p0_recy, 'rx', 'LineWidth',1.5);
    grid on;
    axis equal;
    legend('Original','Reconstructed original');
    title('Shapes in (x,y) space');
    xlabel('x'); ylabel('y');
    
    subplot(1,3,2)
    stem(freq_axis, d0x,'ko-'); hold on
    stem(freq_axis_H, d0Hx,'rx')
    xlabel('Frequency k');
    ylabel('Descriptor d_x');
    title('x descriptors');
    legend('Original','Low frequency');
    grid on;
    
    subplot(1,3,3)
    stem(freq_axis, d0y,'ko-'); hold on
    stem(freq_axis_H, d0Hy,'rx')
    xlabel('Frequency k');
    ylabel('Descriptor d_y');
    title('y descriptors');
    legend('Original','Low frequency');
    grid on;

%% Figure 2: Magnitude of descriptors
figure(2);
    d0_mags  = sqrt(sum(reshape(d0,  2, []).^2, 1));
    d0H_mags = sqrt(sum(reshape(d0H, 2, []).^2, 1));
    
    stem(freq_axis, d0_mags,'o'); hold on
    stem(freq_axis_H, d0H_mags,'x');
    xlabel('Frequency k');
    ylabel('||[d_x,d_y]||');
    title('Descriptor magnitude');
    grid on;

end
