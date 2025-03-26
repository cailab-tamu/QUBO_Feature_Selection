load('../Data_hESC_EC_day1_5000g.mat')
% Flippling embedding coordinates
sce.struct_cell_embeddings.phate2d(:,2) = -sce.struct_cell_embeddings.phate2d(:,2);
sce.struct_cell_embeddings.phate2d(:,1) = -sce.struct_cell_embeddings.phate2d(:,1);
sce.s = sce.struct_cell_embeddings.phate2d;
sce.c = sce.c_cell_type_tx;


f=scgeatool(sce);
%f = figure;
%scatter(sce.s(:,1), sce.s(:,2),10,grp2idx(sce.c_cell_type_tx),"filled");
%colormap(pkg.i_mycolorlines(2));

f.Position(3)=870;
f.Position(4)=320;

filename = 'mapping_pseudotime.png';

% Save as high-quality PNG using print
resolution = 300; 
% print(f, filename, '-dpng', sprintf('-r%d', resolution));
set(gca, 'Fontsize', 15, 'LineWidth', 1.5)
xlabel('PHATE 1')
ylabel('PHATE 2')
box on