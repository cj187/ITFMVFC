function [U,fail,iteration] = ITFMVFC(data,H,W,options,n_cluster,n_data,n_view,d_view,d_h_view,iter_max)
fail = 0;
alpha = options(1);
eta = options(2);
rho = options(3);

select_k= randsample(1:n_view-1,1);
[~,P] = kmeans(data{select_k}',n_cluster,'Replicates', 3);
dist_center_data = pdist2(P,data{select_k}');
dist_center_data(dist_center_data==0) = 1e-6;
U = 1./((dist_center_data.^2).*(ones(n_cluster,1)*sum(dist_center_data.^(-2))));

DW = cell(n_cluster,1);
for i = 1:n_cluster
    % [row_idx, col_idx, ~] = find(sum_W);
    % num_nonzero = length(row_idx);
    % varpi_values = 1e-3 * randn(num_nonzero, 1);
    % 
    % DW{i} = sparse(row_idx, col_idx, varpi_values,n_data,n_data);

    DW{i} = sparse(n_data,n_data);
end

% iter_max = 15;
min_impro = 1e-4;
obj = zeros(iter_max,1);
for iter = 1:iter_max
    U_old = U;
    dist_center_data = cell(n_view,1);
    dist_center_H = cell(n_view,1);
    M = zeros(n_view,1);
    T = zeros(n_view,1);

    for k = 1:n_view-1
        mf = U.^2;
        q = (mf*data{k}')./(sum(mf,2)*ones(1,d_view(k)));
        p = (mf*H{k})./(sum(mf,2)*ones(1,d_h_view));
        dist_center_data{k} = pdist2(q,data{k}');
        dist_center_data{k}(dist_center_data{k}==0) = 1e-6;
        dist_center_H{k} = pdist2(p,H{k});
        dist_center_H{k}(dist_center_H{k}==0) = 1e-6;

        [j_idx, z_idx, w_val] = find(W{k});
        L1_dist = sum(abs(U(:, j_idx) - U(:, z_idx)), 1);
        tmp = alpha/2*sum(w_val .* L1_dist');


        M(k,1) = eta*(sum(sum(mf.*(dist_center_data{k}.^2)))+tmp);
        T(k,1) = (1-eta)*(sum(sum(mf.*(dist_center_H{k}.^2)))+tmp);
    end

    k = n_view;
    mf = U.^2;
    q = (mf*data{k}')./(sum(mf,2)*ones(1,d_view(k)));
    p = (mf*H{k})./(sum(mf,2)*ones(1,d_h_view));
    dist_center_data{k} = pdist2(q,data{k}');
    dist_center_data{k}(dist_center_data{k}==0) = 1e-6;
    dist_center_H{k} = pdist2(p,H{k});
    dist_center_H{k}(dist_center_H{k}==0) = 1e-6;
    M(k,1) = eta*(sum(sum(mf.*(dist_center_data{k}.^2))));
    T(k,1) = (1-eta)*(sum(sum(mf.*(dist_center_H{k}.^2))));
    zeta1 = mean(M);
    zeta2 = mean(T);
    phi = exp(M./(-zeta1));
    phi = phi./sum(phi);

    Psi = exp(T./(-zeta2));
    Psi = Psi./sum(Psi);

    obj(iter) = sum(Psi.*T)+sum(phi.*M)+zeta1*sum(phi.*log(phi))+zeta2*sum(Psi.*log(Psi));

    clear M;clear p;
    clear T;
    clear tmp;
    clear q;

    sum_Psi_phi = eta*phi+(1-eta)*Psi;
    sum_W = sparse(n_data,n_data);
    for k = 1:n_view
        sum_W = sum_Psi_phi(k)*W{k}+sum_W;
    end
    

    n_not_zero = sum(sum_W~=0,2);
    nzero_sum_W = (sum_W~=0);
    epst = 1e-6*nzero_sum_W;
    [x,y,w] = find(sum_W);
    bu = zeros(n_cluster,n_data);
    eu = bu;


    clear sum_W;

    DV = cell(n_cluster,1);
    nn = 1:n_data;
    for i = 1:n_cluster
        cu = sparse(nn,nn,U(i,:)');
        DU = cu*nzero_sum_W;
        sst = rho*abs(DU+1/rho*DW{i}-DU'-1/rho*DW{i}')+epst;
        Ksai = sst(nzero_sum_W);
        Ksai = max(1-alpha*w./Ksai,0.5);
        sst = sparse(x,y,Ksai,n_data,n_data);
        sst1 = sparse(x,y,1-Ksai,n_data,n_data);
        DV{i}  = sst.*(DU+ 1/rho*DW{i})+sst1.*(DU'+ 1/rho*DW{i}');
    end

    Dist = zeros(n_cluster,n_data,n_view);
    parfor k = 1:n_view
        Dist(:,:,k) = eta*phi(k)*dist_center_data{k}.^2+...
            (1-eta)*Psi(k)*dist_center_H{k}.^2;
    end



    parfor i = 1:n_cluster
        bu(i,:) = sum(Dist(i,:,:),3)+ rho/2 *n_not_zero';
        eu(i,:) = sum(DW{i}-rho*DV{i},2);
    end
    clear Dist;
    Kappa = (1+sum(eu./(2*bu),1))./sum(0.5./bu,1);

    parfor i = 1:n_cluster
        U(i,:) = (Kappa-eu(i,:))./(2*bu(i,:));
    end

    parfor i = 1:n_cluster
        cu = sparse(nn,nn,U(i,:)');
        DU = cu*nzero_sum_W;
        DW{i} = DW{i}+rho*(DU-DV{i});
    end

    if iter > 1
        if max(abs(U_old - U),[],'all') < min_impro
        % if abs(obj(iter)-obj(iter-1))/max(obj(iter-1),1)<min_impro
            iteration = iter;
            break;
        end
        iteration = iter;
    end
end
