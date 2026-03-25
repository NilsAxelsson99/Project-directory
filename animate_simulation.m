function animate_simulation(trajectory, time_points, params)
% Animate simulation and save movie
num_steps = size(trajectory,1);

skip = max(1, round(num_steps/200));
%skip = 100*h;
%skip=100;
realTimeScale = 1.0;

N = params.N;
xdim = params.xdim_agent;

% Extract trajectories
pos_traj = zeros(2,N,num_steps);
theta_traj = zeros(N,num_steps);

for k = 1:num_steps
    Xk = trajectory(k,:);
    for i = 1:N
        base = (i-1)*xdim;
        pos_traj(:,i,k) = Xk(base + (1:2));
        theta_traj(i,k) = Xk(base + 3);
    end
end

%% Video setup
video_filename = 'unicycle_simulation.mp4';
v = VideoWriter(video_filename,'MPEG-4');
v.FrameRate = 30;
open(v);

%% Figure setup
figure('Color','w'); hold on; axis equal; grid on;
xlabel('x'); ylabel('y');
title(sprintf('Formation (H=%d)', params.H));

x_all = squeeze(pos_traj(1,:,:));
y_all = squeeze(pos_traj(2,:,:));
padding = 1.0;
xlim([min(x_all(:))-padding, max(x_all(:))+padding]);
ylim([min(y_all(:))-padding, max(y_all(:))+padding]);

trailLines   = gobjects(N,1);
agentMarkers = gobjects(N,1);
headingLines = gobjects(N,1);
agentLabels  = gobjects(N,1);

for i = 1:N
    trailLines(i)   = plot(nan,nan,'-','Color',[0.6 0.8 1],'LineWidth',1.2);
    agentMarkers(i) = plot(nan,nan,'ko','MarkerFaceColor','b','MarkerSize',8);
    headingLines(i) = plot(nan,nan,'k-','LineWidth',1.0);

    % Agent number label
    agentLabels(i) = text(nan, nan, sprintf('%d', i), ...
        'FontSize', 10, ...
        'FontWeight', 'bold', ...
        'Color', 'k', ...
        'HorizontalAlignment', 'left', ...
        'VerticalAlignment', 'bottom');
end

outline_handle = plot(nan,nan,'r--','LineWidth',1.2);

centroid_true_handle = plot(nan,nan,'rx','MarkerSize',14,'LineWidth',2);
centroid_ref_handle  = plot(nan,nan,'go','MarkerSize',10,'LineWidth',2);

%% Animation loop
for k = 1:skip:num_steps
    for i = 1:N
        xs = squeeze(pos_traj(1,i,1:k));
        ys = squeeze(pos_traj(2,i,1:k));
        set(trailLines(i),'XData',xs,'YData',ys);
        set(agentMarkers(i),'XData',xs(end),'YData',ys(end));

        % Slight offset so number does not overlap marker
    offset = 0.15;
    set(agentLabels(i), ...
    'Position', [xs(end)+offset, ys(end)+offset, 0]);

        th = theta_traj(i,k);
        L = 0.4;
        set(headingLines(i), ...
            'XData',[xs(end), xs(end)+L*cos(th)], ...
            'YData',[ys(end), ys(end)+L*sin(th)]);
    end

    Xkpos = squeeze(pos_traj(:,:,k));
    set(outline_handle, ...
        'XData',[Xkpos(1,:), Xkpos(1,1)], ...
        'YData',[Xkpos(2,:), Xkpos(2,1)]);

    centroid_true = mean(Xkpos,2);
    centroid_ref  = params.centroid_ref(time_points(k));

    set(centroid_true_handle,'XData',centroid_true(1),'YData',centroid_true(2));
    set(centroid_ref_handle, 'XData',centroid_ref(1), 'YData',centroid_ref(2));

    drawnow limitrate;
    axis equal

    writeVideo(v, getframe(gcf));
    %pause(skip*params.h/realTimeScale);
    %pause
    %dt = time_points(min(k+skip,end)) - time_points(k);
    %pause(dt / realTimeScale);
end

close(v);
disp(['Saved movie to: ', video_filename]);

end
