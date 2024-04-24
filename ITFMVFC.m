function [U,fail,iteration] = ITFMVFC(data,H,W,options,n_cluster,n_data,n_view,d_view,eps)
fail = 0;
alpha = options(1);
zeta1 = options(2);
zeta2 = options(3);
eta = options(4);
rho = options(5);
U = zeros(n_cluster, n_data);
for k = 1:n_view
    [~,P0{k}] = kmeans(data{k}',n_cluster);
    dist_center_data{k} = pdist2(P0{k},data{k}');
    dist_center_data{k}(dist_center_data{k}==0) = 1e-6;
end
U = randn(n_cluster,n_data);
iter_max = 150;
min_impro = eps;

for iter = 1:iter_max
    U_old = U;
    q = cell(n_view,1);
    p = q;
    dist_center_data = q;
    dist_center_H = q;
    M = zeros(n_view,1);
    T = zeros(n_view,1);

    for k = 1:n_view
        mf = U.^2;
        q{k} = (mf*data{k}')./(sum(mf,2)*ones(1,d_view(k)));
        p{k} = (mf*H{k})./(sum(mf,2)*ones(1,n_cluster));
        dist_center_data{k} = pdist2(q{k},data{k}');
        dist_center_data{k}(dist_center_data{k}==0) = 1e-6;
        dist_center_H{k} = pdist2(p{k},H{k});
        dist_center_H{k}(dist_center_H{k}==0) = 1e-6;

        tmp = alpha/2*(sum(sum(W{k}.*(squareform(pdist(U','minkowski',1))))));
        M(k,1) = eta*(sum(sum(mf.*(dist_center_data{k}.^2)))+tmp);
        T(k,1) = (1-eta)*(sum(sum(mf.*(dist_center_H{k}.^2)))+tmp);
    end
    phi = exp(M./(-zeta1));
    phi = phi./sum(phi);

    Psi = exp(T./(-zeta2));
    Psi = Psi./sum(Psi);

    if(any(isnan(phi)) || any(isnan(Psi)))
        fail = 1;
        iteration = 0;
        break;
    end

    sum_Psi_phi = eta*phi+(1-eta)*Psi;
    sum_W = zeros(n_data,n_data);
    for k = 1:n_view
        sum_W = sum_Psi_phi(k)*W{k}+sum_W;
    end

    n_not_zero = sum(sum_W~=0,2);
    nzero_sum_W = (sum_W~=0);
    epst = 1e-4*nzero_sum_W;
    [x,y,w] = find(sum_W);
    au = zeros(n_cluster,n_data);
    bu = au;

    DW = cell(n_cluster,1);
    for i = 1:n_cluster
        DW{i} = sum_W;
    end

    DV = DW;
    nn = 1:n_data;
    for i = 1:n_cluster
        cu = sparse(nn,nn,U(i,:)');
        DU{i} = cu*nzero_sum_W;
        sst = rho*abs(DU{i}+1/rho*DW{i}-DU{i}'-1/rho*DW{i}'+epst);
        Ksai = sst(nzero_sum_W);
        Ksai = max(1-alpha*w./Ksai,0.5);
        sst = sparse(x,y,Ksai);
        sst1 = sparse(x,y,1-Ksai);
        DV{i}  = sst.*(DU{i}+ 1/rho*DW{i})+sst1.*(DU{i}'+ 1/rho*DW{i}');
        DW{i} = DW{i}+rho*(DU{i}-DV{i});
    end

    for k = 1:n_view
        Dist(:,:,k) = eta*phi(k)*dist_center_data{k}.^2+...
            (1-eta)*Psi(k)*dist_center_H{k}.^2;
    end

    for i = 1:n_cluster
        au(i,:) = sum(Dist(i,:,:),3)+ rho/2 *n_not_zero';
        bu(i,:) = sum(DW{i}-rho*DV{i},2);
    end
    Kappa = (1+sum(bu./(2*au),1))./sum(0.5./au,1);

    for i = 1:n_cluster
        U(i,:) = (Kappa-bu(i,:))./(2*au(i,:));
    end

    if iter > 1
        if max(abs(U_old - U)) < min_impro
            iteration = iter;
            break;
        end
        iteration = iter;
    end
end
