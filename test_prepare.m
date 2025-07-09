clear;clc;
close all;

load('./SYN_1_Ori.mat');


[uniqueLabels, ~, labelIndices] = unique(truelabel);
samplesPerCluster = accumarray(labelIndices, 1);

rand_i = 3;
axi_max = 14;

rng(rand_i);
numClusters = length(samplesPerCluster);

centers(1, :) = rand(1, 2) * axi_max;
minDist = (axi_max/numClusters)*0.8;

% 迭代选择后续的聚类中心
for i = 2:numClusters
    tooClose = true;
    while tooClose
        candidateCenter = rand(1, 2) * axi_max;
        dists = pdist2(centers(1:i-1, :),candidateCenter);
        % dists = sqrt(sum((centers(1:i-1, :) - candidateCenter).^2, 2));
        if min(dists(:))>minDist
            tooClose = false;
            centers(i, :) = candidateCenter;
        else
            fprintf('too close\n');
        end
    end
end

tmp_data = [];
% 随机生成每个簇的数据
for i = 1:numClusters
    % 对于每个簇，随机选择一个均值和协方差
    % meanVal = rand(1, 2) * 40; % 随机均值在0到20之间
    meanVal = centers(i, :);
    covVal = diag(rand(1, 2) * 2); % 随机方差在0到2之间
    
    % 使用mvnrnd生成正态分布数据
    clusterData = mvnrnd(meanVal, covVal, samplesPerCluster(i));
    tmp_data = [tmp_data; clusterData];
end
data_view2 = tmp_data;

data_view3 = [sin(data_view2(:,1)), log1p(abs(data_view2(:,2)))];
