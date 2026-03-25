function [x0_agent, p0] = init_agents(params)

N = params.N;
spread = 5; % 5; 

p0 = randn(2,N)*spread;
theta0 = 2*pi*rand(1,N) - pi;

x0_agent = zeros(N*params.xdim_agent,1);
    for i = 1:N
        base = (i-1)*params.xdim_agent;
        x0_agent(base+1:base+3) = [p0(:,i); theta0(i)];
    end
end
