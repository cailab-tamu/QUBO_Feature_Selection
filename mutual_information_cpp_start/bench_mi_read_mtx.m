% Clear workspace
clc;
clearvars;

% Specify the file name
filename = 'sparse_matrix.mtx';

% Open the file
fid = fopen(filename, 'r');

% Read the header lines (you can skip comments starting with '%' or additional metadata)
while true
    line = fgetl(fid);
    if line(1) ~= '%'
        break;
    end
end

% Read the data (assuming the file has 3 columns: row, column, value)
data = textscan(fid, '%f %f %f');

% Close the file
fclose(fid);

% Convert data to a sparse matrix
rows = data{1}; % Row indices
cols = data{2}; % Column indices
vals = data{3}; % Values

% Determine the size of the sparse matrix
numRows = max(rows);
numCols = max(cols);

% Create the sparse matrix
sparseMatrix = sparse(rows, cols, vals, numRows, numCols);

tic;
MI_mat = MI_block_construction(sparseMatrix);
time_end = toc;
fprintf("MI computation time with pool-start: %3.4f \n",time_end);

tic;
MI_mat = MI_block_construction(sparseMatrix);
time_end = toc;
fprintf("MI computation time : %3.4f \n",time_end);