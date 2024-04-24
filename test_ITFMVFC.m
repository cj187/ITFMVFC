clc;
clear;

root_path = 'C:/Users/Documents/MATLAB/ITFMVFC/';


%% load data
matrix_data_1 = {'1,3sources','2,Caltech101-7','3,Caltech_2','4,motion','5,WebKB','6,prokaryotic','7,100Leaves','8,20newsgroups','9,uci-digit',...
    '10,MSRC','11,Caltech101-20'};
matrix_data_2 = {'1,SYN_1','2,SYN_2'};

%% Parameters setting
datatype = 2;%1 real-world, 2 synthetic
n_index = 12;%评价指标的个数，包括迭代次数和运行时间

alpha_candidates = 0.01; %[0.001,0.01,0.1,1,10,100,1000]
zeta_candidates = 1e3; %[1e2,1e3,1e4,1e5,1e6]
eta_candidates = 0.6; %[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3,0.2,0.1,0]
rho_candidates = 2; %[0.01,0.1,1,2,4,10]
eps = 1e-4;

re_times = 20;
algName = 'ITF-MVFC';


if datatype==1
    matrix_data = matrix_data_1;
    type_str = 'RealworldData3/';
    Q = 20;
elseif datatype == 2
    matrix_data = matrix_data_2;
    type_str = 'synthetic/';
    Q = 20;
end

for ii = 1:1
    for data_i = 1:2
        dataname = matrix_data{data_i};
        dataname = strsplit(dataname,',');
        dataname = dataname{2};
        load(['./synthetic/',dataname,'.mat']);

        n_view = length(data);
        n_data = length(truelabel);
        n_cluster = max(truelabel);
        d_h_view = n_cluster;
        d_h_tview = n_cluster;
        F = cell(n_view,1);
        H = cell(n_view,1);
        W = topology;
        %% data的标准化
        for k = 1:n_view
            data{k} = normalize_data(data{k}');
            data{k} = data{k}';
            d_view(k) = size(data{k},1);
            F{k} = rand(d_view(k),d_h_view);
            H{k} = rand(n_data,d_h_tview);
        end
        %% NMF求隐形物理特征空间B,隐形拓扑特征空间H和S
        B = rand(d_h_view,n_data);
        myeps = 1e-8;
        max_iter = 200;
        q = rand(n_view,1);
        q = q./sum(q);

        S = rand(d_h_tview,n_data);
        q2 = rand(n_view,1);
        q2 = q2./sum(q2);
        for iter = 1:max_iter
            numerator_B = zeros(d_h_view,n_data);
            denominator_B = zeros(d_h_view,n_data);
            numerator_S = zeros(d_h_tview,n_data);
            denominator_S = numerator_S;
            for k = 1:n_view
                F{k} = F{k}.*((data{k}*B')./max(F{k}*(B*B'),myeps));
                numerator_B = numerator_B+F{k}'*data{k}*q(k);
                denominator_B = denominator_B+(F{k}'*F{k})*B*q(k);

                H{k} = H{k}.*((W{k}*S')./max(H{k}*(S*S'),myeps));
                numerator_S = numerator_S+H{k}'*W{k}*q2(k);
                denominator_S = denominator_S+(H{k}'*H{k})*S*q2(k);
            end
            B = B.*(numerator_B./max(denominator_B,myeps));
            S = S.*(numerator_S./max(denominator_S,myeps));

            for k = 1:n_view
                tmp(k,1) = norm(data{k} - F{k} * B)*(-1/1e3);
                tmp2(k,1) = norm(W{k} - H{k} * S)*(-1/1e3);
            end
            q = tmp./sum(tmp);
            q2 = tmp2./sum(tmp2);
        end
        data{n_view+1} = B;
        d_view(n_view+1) = size(B,1);
        H{n_view+1} = S';
        n_view = n_view+1;
        for k = 1:n_view
            d_view(k) = size(data{k},1);
        end
        W{n_view} = zeros(n_data,n_data);


        for model_i = 5:5
            fprintf('Begin: Data %d model %d ii %d\n', data_i,model_i,ii);
            model = model_i;

            for alpha_i = 1:length(alpha_candidates)
                alpha = alpha_candidates(alpha_i);
                for zeta_i = 1:length(zeta_candidates)
                    zeta2 = zeta_candidates(zeta_i);
                    zeta1 = zeta_candidates(zeta_i);
                    for eta_i = 1:length(eta_candidates)
                        eta = eta_candidates(eta_i);
                        for rho_i = 1:length(rho_candidates)
                            rho = rho_candidates(rho_i);
                            options = [alpha,zeta1,zeta2,eta,rho];
                            rerun_i = 1;
                            while rerun_i <= re_times
                                [U,fail,tmp_iter] = DVHVC(data,H,W,options,n_cluster,n_data,n_view,d_view,eps);
                            end
                        end
                    end
                end
            end
        end
    end
