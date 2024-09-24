function MI_mat_block = MI_block(data_block, data_block2, mode, comp_ii)
    % MI_construction computes the mutual information from data
    % INPUT:
    % data =====> Contains X count matrix and y target as follows
    %             data -> (cells + target) by genes 
    % mode =====> mode 1 is upper triangular, mode 2 is full
    % comp_ii ==> Compute diagonal terms? true/false
    % USAGE:
    % R0 = MI_construction(data);
    % save('R0.mat', 'R0', '-v7.3');

    if nargin < 2 || isempty(data_block2); data_block2 = data_block; end
    if nargin < 3 || isempty(mode); mode = 1; end
    if nargin < 4 || isempty(comp_ii); comp_ii = false; end 

    offset = 1;
    % remove offset
    if comp_ii 
        offset = 0;
    end
    % Data is count matrix with genes in rows and cells in columns
    %data = sparse(data);
    %data2 = sparse(data2);
    % transpose for efficient access (vectorization) across observations
    %data = data';
    %data2 = data2';

    nobs = size(data_block,1);
    nobs2 = size(data_block2, 1);
    if nobs ~= nobs2
        error("Not same number of observations in input matrices");
    end
    ngene = size(data_block,2);
    ngene2 = size(data_block2, 2);
    MI_mat_block = zeros(ngene, ngene2);
    % MI across data's rows (Computing upper triangular) in parallel 
    if mode == 1
        parfor ig = 1:ngene
            tmp = zeros(ngene2, 1);
            for jg = ig + offset:ngene2 % Skip the ig-ig term
                % Computing pair MI with binning distribution
                tmp(jg) = BinPairMI( full(data_block(1:nobs, ig)), ...
                                     full(data_block2(1:nobs2, jg)) );
            end
            MI_mat_block(ig, :) = tmp;
        end
    else
        parfor ig = 1:ngene
            tmp = zeros(ngene2,1);
            for jg = 1:ngene2
                % Computing pair MI with binning distribution
                tmp(jg) = BinPairMI( full(data_block(1:nobs, ig)), ...
                                     full(data_block2(1:nobs2, jg)) );
            end
            MI_mat_block(ig,:) = tmp;
        end
    end
end