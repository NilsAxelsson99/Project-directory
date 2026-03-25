function Neighbors = build_ring_graph(N)

Neighbors = cell(N,1);
    for i = 1:N
        Neighbors{i} = mod([i-2, i], N) + 1;
    end
end
