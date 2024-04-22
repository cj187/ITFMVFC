function [ ids ] = new_spectral_clustering(W, numClusters)

D = diag(1./sqrt(sum(W, 2)+ eps));
W = D * W * D;
[U, s, V] = svd(W);
V = U(:, 1 : numClusters);
% if ~isreal(V)
%     V = V ./ (sum((V .* conj(V)).^2)).^0.25;
% end
% V = normr(V);
if isreal(V)
    V = normr(V);
end

ids = litekmeans(V, numClusters, 'MaxIter',100, 'Replicates',200);

end
