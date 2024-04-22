function data_ori = dataNormalization(data_ori,model)
switch model
    case 1
        data_ori=data_normalize(data_ori,'range');
        data_ori=normr(data_ori);
    case 2
        data_ori=data_normalize(data_ori,'var');
        data_ori=normr(data_ori);
    case 3
        data_ori=data_normalize(data_ori,'xy');
        data_ori=normr(data_ori);
    case 4
        data_ori=normr(data_ori);
    case 5
        data_ori=data_normalize(data_ori,'range');
    case 6
        data_ori=data_normalize(data_ori,'var');
    case 7
        data_ori=data_normalize(data_ori,'xy');
    case 8
    otherwise
        warning('error normalization.');
end