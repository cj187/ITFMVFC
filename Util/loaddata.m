function [data, W ,ground_truth] = loaddata(dataname,datatype,root_path)

switch datatype
    case 1
        load([root_path,'RealworldData2/',dataname,'.mat']);
        data_ori = data;
        W_ori = topology;
        ground_truth_ori = truelabel;
%         ground_truth_ori = load([root_path,'data/synthetic/label/label_',dataname,'.txt']);
%         ground_truth_ori = 9*ones(size(data_ori,1),1);
    case 2
        load([root_path,'ArtificialDataset/',dataname,'.mat']);
        data_ori = data;
        W_ori = topology;
        ground_truth_ori = truelabel;
    case 3
        data_ori = load([root_path,'data/gene/clean_data/data_',dataname,'.txt']);
        ground_truth_ori = load([root_path,'data/gene/label/label_',dataname,'.txt']);
    case 4
        data_ori = load([root_path,'data/GEPFCMdata/data_',dataname,'.txt']);
        ground_truth_ori = load([root_path,'data/GEPFCMdata/label_',dataname,'.txt']);
    otherwise
        error('Please check data type\n');
end

data = data_ori;
W = W_ori;
ground_truth = ground_truth_ori;
% [data,remain_dx,~] = unique(data_ori,'rows','stable');
% ground_truth = ground_truth_ori(remain_dx);
% if size(data,1)~=size(data_ori,1)
%     fprintf('There are equal rows, original number of data is %d\n', size(data_ori,1));
% end