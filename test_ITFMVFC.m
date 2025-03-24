clc;
clear;

load('./Datasets/SYN_1.mat');

n_view = length(data);
n_data = length(truelabel);
n_cluster = max(truelabel);

W = cell(n_view,1);
for k = 1:n_view
    tmp_data = data{k}';
    [j_idx, z_idx, ~] = find(topology{k});
    distances_sq = sum((tmp_data(j_idx,:) - tmp_data(z_idx,:)).^2, 2);
    tmp = exp(-distances_sq / 0.5);
    W{k} = sparse(j_idx, z_idx, tmp, size(topology{k},1), size(topology{k},2));
end


d_view=zeros(n_view,1);
for k = 1:n_view
    d_view(k) = size(data{k},1);
end

d_h_view = min(n_cluster,min(d_view));
d_h_tview = d_h_view;

max_iter = 50;
[F, B,iter] = multi_view_nmf(data, d_h_view, max_iter);
[H, S,iter2] = multi_view_nmf(W, d_h_tview, max_iter);

data{n_view+1} = B;
d_view(n_view+1) = size(B,1);
H{n_view+1} = S';

n_view = n_view+1;
d_view(n_view) = d_h_view;
W{n_view} = sparse(n_data,n_data);


rho = 1.0;
alpha=2^(3);
eta=0.5;
zeta1 = 1000;
zeta2 = 1000;
iter_max = 100;
options = [alpha, eta, rho,zeta1,zeta2];

[U,fail,tmp_iter] = ITFMVFC(data,H,W,options,n_cluster,...
n_data,n_view,d_view,d_h_view,iter_max);





