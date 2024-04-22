function w =  comput_hw(S,k)
n_data = size(S,1);
dist_data = pdist(S);
dist_data = dist_data.^2;
dist_data = squareform(dist_data);
w1 = zeros(n_data,n_data);
w = zeros(n_data,n_data);
for i = 1:n_data
    tmp_dist = dist_data(i,:);
    [~,orddist] = sort(tmp_dist);
    w1(i,orddist(2:k+1)) = 1;
    tmp_w = exp(-tmp_dist./0.5);
    w(i,:) = w1(i,:).*tmp_dist;
end
