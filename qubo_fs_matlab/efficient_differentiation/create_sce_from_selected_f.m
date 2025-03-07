glnot = Tml{4}.selectedGenes;
gl = Tqubo.selectedGenes;

gln = intersect(gl,glnot,'stable');

gidx = ismember(sce.g, gln);
g = sce.g(gidx);
X = sce.X(gidx,:);

sce2 = SingleCellExperiment(X, g);      % make SCE class
sce2.c_cell_id = sce.c_cell_id;
sce2.c_batch_id = sce.c_batch_id;
sce2.c_cell_type_tx = sce.c_cell_type_tx;
sce2 = leiden_clustering_ann(sce2, 0.5);
scgeatool(sce2)