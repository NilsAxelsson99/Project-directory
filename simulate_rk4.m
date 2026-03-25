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

