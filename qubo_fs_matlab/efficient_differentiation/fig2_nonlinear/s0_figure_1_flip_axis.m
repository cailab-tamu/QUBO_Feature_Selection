load('../D1_1f_G2M_5000g_reMonocle.mat')
% Flippling embedding coordinates
sce.struct_cell_embeddings.phate2d(:,2) = -sce.struct_cell_embeddings.phate2d(:,2);
sce.struct_cell_embeddings.phate2d(:,1) = -sce.struct_cell_embeddings.phate2d(:,1);
sce.s = sce.struct_cell_embeddings.phate2d;
f=scgeatool(sce);
f.Position(3)=450;
f.Position(4)=250;
box on