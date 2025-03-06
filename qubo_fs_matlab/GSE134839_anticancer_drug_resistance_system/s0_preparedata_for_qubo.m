load('GSE134839_cleandata.mat')
idx = find(strcmp('manual_pseudotime',sce.list_cell_attributes(1:2:end)));
t = sce.list_cell_attributes{idx+1};



