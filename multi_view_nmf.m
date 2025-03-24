function [F, B, iter] = multi_view_nmf(X, r, max_iter)
    K = length(X); % 视图数
    n = size(X{1}, 2); % 样本数
    F = cell(K, 1);  % 存储每个视图的F
    d = zeros(K, 1); % 每个视图的特征维度

    % 初始化 F 和 B
    for k = 1:K
        d(k) = size(X{k}, 1);
        F{k} = rand(d(k), r);  % 初始化每个视图的基矩阵
    end
    B = rand(r, n);  % 共享表示矩阵

    for iter = 1:max_iter
        prev_B = B; % 用于跳出循环的判断
        % **更新 F{k}**
        for k = 1:K
            F{k} = F{k} .* ((X{k} * B') ./ max(F{k} * (B * B'), eps));
        end
        
        % **更新 B**
        numerator = zeros(r, n);
        denominator = zeros(r, n);
        for k = 1:K
            numerator = numerator + F{k}' * X{k};
            denominator = denominator + F{k}' * F{k} * B;
        end
        hasNann = any(isnan(numerator), 'all');
        hasNand = any(isnan(denominator), 'all');
        if hasNand | hasNann
            t=12;
        end
        B = B .* (numerator ./ max(denominator, eps));

        hasNanB = any(isnan(B), 'all');
        if hasNanB
            t =12;
        end


        if norm(B - prev_B, 'fro') / norm(prev_B, 'fro') < 1e-3
            break;
        end
        % prev_B = B;
    end

    fprintf('iter=%d \t %.5f\n',iter, norm(B - prev_B, 'fro') / norm(prev_B, 'fro'));
end
