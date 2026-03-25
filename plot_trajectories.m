function plot_trajectories(trajectory, params)
% Plot agent trajectories in the plane

N = params.N;
xdim = params.xdim_agent;
num_steps = size(trajectory,1);

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

figure;
hold on; grid on; axis equal;
xlabel('x'); ylabel('y');
title('Agent trajectories');

    for i = 1:N
        plot(squeeze(pos_traj(1,i,:)), ...
             squeeze(pos_traj(2,i,:)), ...
             'LineWidth',1.2);
    end

end
