function saveResult(path_process,path_analysis,info,cell_index)

for i_index = 1:length(cell_index)
    tmp_index = cell_index{i_index};
    t = table(sort(roundn(tmp_index,-4))',info,datetime('now'));
    writetable(t,path_process,Sheet=i_index,WriteMode='append',WriteRowNames=false,WriteVariableNames=false);

    t = table(roundn(max(tmp_index),-4),roundn(mean(tmp_index),-4),roundn(min(tmp_index),-4),roundn(std(tmp_index),-4),datetime('now'));
    writetable(t,path_analysis,Sheet=i_index,WriteMode='append',WriteRowNames=false,WriteVariableNames=false);
end