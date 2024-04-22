clc;
clear;

root_path = 'C:/Users/18738/Documents/MATLAB/DVHVC_MVFC/';
% root_path = 'D:/Matlab/exp_cj/';

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

                            %1:ACC; 2:NMI; 3:ARI; 4:F; 5:P; 6:R; 7:Purity; 8:CHI; 9:DBI; 10:SCI; 11:iter; 12:time
                            index_name = ["ACC";"NMI";"ARI";"F";"P";"R";"Purity";"CH";"DB";"SC";"iter";"time"];
                            cell_index = cell(n_index,1);

                            rerun_i = 1;
                            err_times = 0;
                            err_times1 = 0;
                            isSave = 0;

                            while rerun_i <= re_times
                                tic;
                                [U,fail,tmp_iter] = DVHVC(data,H,W,options,n_cluster,n_data,n_view,d_view,eps);
                                if(fail == 1)
                                    err_times1 = err_times1+1;
                                    if err_times1> 5
                                        fprintf('Fail Error\n');
                                        break;
                                    end
                                else
                                    err_times1 = 0;
                                    cell_index{12}(rerun_i,1) = toc;
                                    cell_index{11}(rerun_i,1) = tmp_iter;
                                    [~,label_all] = max(U);
                                    idx = bestMap(truelabel,label_all);

                                    cell_index{1}(rerun_i,1) = sum(idx==truelabel)/length(label_all);
                                    cell_index{2}(rerun_i,1) = compute_nmi(truelabel,idx); % NMI
                                    cell_index{3}(rerun_i,1) = compute_rand_index(truelabel,idx); % ARI
                                    [cell_index{4}(rerun_i,1),cell_index{5}(rerun_i,1),cell_index{6}(rerun_i,1)] = compute_f(truelabel,idx); % F1, precision, recall
                                    cell_index{7}(rerun_i,1) = compute_purity(truelabel,idx); %PUR

                                    tmp_index = cell(3,1);
                                    for k = 1:n_view
                                        tmp_index{1}(k,1) = evalclusters(data{k}',idx,'CalinskiHarabasz').CriterionValues;
                                        tmp_index{2}(k,1) = evalclusters(data{k}',idx,'DaviesBouldin').CriterionValues;
                                        tmp_index{3}(k,1) = evalclusters(data{k}',idx,'silhouette').CriterionValues;
                                    end

                                    for i_index = 1:3
                                        cell_index{i_index+7}(rerun_i,1) = mean(tmp_index{i_index});
                                    end

                                    if ~isnan(cell_index{8}(rerun_i,1))&&~isnan(cell_index{9}(rerun_i,1))&&~isnan(cell_index{10}(rerun_i,1))
                                        rerun_i = rerun_i+1;
                                        isSave = 1;
                                        err_times = 0;
                                    else
                                        err_times = err_times+1;
                                        if err_times > 5
                                            fprintf('No proper initialization! rerun_i=%d\n',rerun_i);
                                            for i_index = 1:n_index
                                                cell_index{i_index}(rerun_i,:) = [];
                                            end
                                            break;
                                        end

                                    end


                                end

                            end

                            if isSave==1
                                %% 分析指标
                                max_index = zeros(n_index,1);mean_index = zeros(n_index,1);min_index = zeros(n_index,1);
                                for i_index = 1:n_index
                                    tmp_index = cell_index{i_index};
                                    max_index(i_index) = max(tmp_index);
                                    mean_index(i_index) = mean(tmp_index);
                                    min_index(i_index) = min(tmp_index);
                                end
                                max_index = [index_name,num2str(max_index)];
                                mean_index = [index_name,num2str(mean_index)];
                                min_index = [index_name,num2str(min_index)];

                                %% 保存指标
                                folder2 = strcat('./Result/',algName,'/',type_str,num2str(data_i),'_',dataname,'/');
                                if ~exist(folder2,'dir')
                                    mkdir(folder2);
                                end
                                path_process = strcat(folder2,'process_result_',algName,'_',dataname,'_', num2str(model),'.xlsx');
                                info = {'alpha',alpha,'zeta1',zeta1,'zeta2',zeta2,'eta',eta,'rho',rho};

                                path_analysis = strcat(folder2,'analysis_result_',algName,'_',dataname,'_', num2str(model),'.xlsx');

                                saveResult(path_process,path_analysis,info,cell_index);
                                fprintf('Save Result Success!\n');
                            end
                            fprintf('End: Data %d model %d ii %d\n', data_i,model_i,ii);
                        end
                    end
                end
            end
        end
    end
end

