function [G,n_same,n_diff] = generate_edges_prob(N, labels, num_edges, p_intra, p_inter, isRetry)
    %rng('shuffle');
    s = zeros(num_edges,1);
t = zeros(num_edges,1);
added_edges = 0;
n_same = 0;
n_diff = 0;

max_trials = 1e7;
trials = 0;
edge_check = sparse(N, N);

if isRetry
while added_edges < num_edges && trials < max_trials
    i = randi(N);
    j = randi(N);
    if i == j|| edge_check(i, j) > 0 || edge_check(j, i) > 0
        trials = trials + 1;
        continue;
    end
    is_same = labels(i) == labels(j);
    prob = is_same * p_intra + (~is_same) * p_inter;
    if rand < prob
        added_edges = added_edges + 1;
        s(added_edges) = i;
        t(added_edges) = j;
        edge_check(i, j) = 1;
        edge_check(j, i) = 1;  % 无向图对称
        if is_same
            n_same =n_same+1;
        else
            n_diff = n_diff+1;
        end
    else
        trials = trials + 1;
    end
    
end
    if trials == max_trials
        warning('Maximum sampling attempts reached. The required number of edges may not have been fully satisfied');
    end
else
    while added_edges < num_edges
        i = randi(N);
        j = randi(N);
        if i == j
            continue;
        end
        is_same = labels(i) == labels(j);
        prob = is_same * p_intra + (~is_same) * p_inter;
        if rand < prob
            added_edges = added_edges + 1;
            s(added_edges) = i;
            t(added_edges) = j;
            edge_check(i, j) = 1;
            edge_check(j, i) = 1;  % 无向图对称
            if is_same
                n_same =n_same+1;
            else
                n_diff = n_diff+1;
            end
        end
        
    end
end

G = graph(s, t, rand(num_edges, 1), N);
end
