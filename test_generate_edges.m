clc;clear;

load('./Datasets/SYN_1.mat');

n_view = length(data);
n_data = length(truelabel);

for k = 1:n_view
    data{k} = data{k}';
end

isRetry = 0;%Whether to continue sampling iteratively to approach the target number of edges as closely as possible.
sparsity_factor = 10.0/n_data;
p_inter = 0.3;
p_intra = 0.8;
num_edges = floor((n_data * (n_data - 1) / 2) * sparsity_factor);
n_same = zeros(n_view,1);
n_diff = zeros(n_view,1);
W = cell(n_view,1);
topology = cell(n_view,1);
for k = 1:n_view         
    [G,n_same(k),n_diff(k)] = generate_edges_prob(n_data, truelabel, num_edges, p_intra, p_inter, isRetry);
    
    topology{k} = adjacency(G);
    [j_idx, z_idx, ~] = find(topology{k});
    distances_sq = sum((data{k}(j_idx,:) - data{k}(z_idx,:)).^2, 2);
    tmp = exp(-distances_sq / 0.5);
    W{k} = sparse(j_idx, z_idx, tmp, size(topology{k},1), size(topology{k},2));
end
