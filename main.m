%% add to the path of MRST toolbox
startup
clc;
%% Define the geometry
G = cartGrid([30 30],[2 2]);
G = computeGeometry(G);
% Remove the circular holes from the grid
r1 = sum(bsxfun(@minus,G.cells.centroids,[0.25 0.25]).^2,2);
r2 = sum(bsxfun(@minus,G.cells.centroids,[0.375 1]).^2,2);
r3 = sum(bsxfun(@minus,G.cells.centroids,[0.25 1.75]).^2,2);
r4 = sum(bsxfun(@minus,G.cells.centroids,[1 1.625]).^2,2);
r5 = sum(bsxfun(@minus,G.cells.centroids,[1.75 1.75]).^2,2);
r6 = sum(bsxfun(@minus,G.cells.centroids,[1.625 1]).^2,2);
r7 = sum(bsxfun(@minus,G.cells.centroids,[1.75 .25]).^2,2);
G = extractSubgrid(G, (r1>0.015625) & (r2>0.015625) & (r3>0.015625) & ...
    (r4>0.015625) & (r5>0.015625) & (r6>0.015625) & (r7>0.015625));
% Remove the diamond shape from the grid [diamond]
center = [1, 0.646];  % Center coordinates of the square
side_length= 0.5;     % Size of the square (side length)
angle = 45;           % Rotation angle in degrees


% Compute the vertices of the rotated square
half_length = side_length / 2;
rotationMatrix = [cosd(angle), -sind(angle); sind(angle), cosd(angle)];
squareVertices = center + half_length * [-1, -1; -1, 1; 1, 1; 1, -1] * rotationMatrix;

excludecells = false(G.cells.num, 1);

% Iterate through grid cells and exclude those inside the square
for cellIndex = 1:G.cells.num
    % Get the cell centroid
    centroid = G.cells.centroids(cellIndex, 1:2);
    
    % Check if the centroid is inside the rotated square
    if inpolygon(centroid(1), centroid(2), squareVertices(:, 1), squareVertices(:, 2))
        excludecells(cellIndex) = true;
    end
end

% Remove the cells that intersect with the rotated square
G = removeCells(G, excludecells);
% Visualize the final grid
figure;
plotGrid(G);
title('Final Grid'); 
%% Boundary Cells
% Number of cells in x direction
Nx = G.cartDims(1);
cell_size = 2/Nx; % size of each cell

% Find the left boundary cells
left_cells = G.cells.centroids(:, 1) < cell_size;
left_cells_ind = find(left_cells);
% Plot and indicate the left boundary cells
cell_colors = ones(G.cells.num, 1); % Initialize with ones
cell_colors(left_cells_ind) = 0; % Change to zero for left cells
figure;
plotCellData(G, cell_colors);
title('Left Boundary Cells');

% Plot and indicate the right boundary cells
right_cells = G.cells.centroids(:, 1) > 2-cell_size;
right_cells_ind = find(right_cells);
% Plot and indicate the right boundary cells
cell_colors = ones(G.cells.num, 1); % Initialize with ones
cell_colors(right_cells_ind) = 0; % Change to zero for right cells
figure;
plotCellData(G, cell_colors);
title('Right Boundary Cells'); 
%% Newton loop for convergence
nc = G.cells.num;
T = initVariablesADI(ones(nc, 1));
[tol, incr] = deal(1e-4, 1e4);
while norm(incr)>tol
      res = residual(G,T,left_cells_ind,right_cells_ind);
      % Newton update
      incr = -res.jac{1}\res.val;
      % New Guess for next iteration
      T = T + incr;
end
%% Plotting and visualization
figure;
plotCellData(G,T.val);
title('Temperature Field'); 
display(T.val)