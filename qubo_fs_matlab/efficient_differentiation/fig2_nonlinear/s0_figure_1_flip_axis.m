load('../D1_1f_G2M_5000g_Monocle.mat')
% Flippling embedding coordinates
sce.struct_cell_embeddings.phate2d(:,2) = -sce.struct_cell_embeddings.phate2d(:,2);
sce.struct_cell_embeddings.phate2d(:,1) = -sce.struct_cell_embeddings.phate2d(:,1);
sce.s = sce.struct_cell_embeddings.phate2d;
sce.c = sce.c_cell_type_tx;
f=scgeatool(sce);
f.Position(3)=870;
f.Position(4)=300;
box on