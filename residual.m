%% Residual Function
function res = residual(G,T,left_boundary,right_boundary)
    % Discreate Operators
    C = G.faces.neighbors;
    C = C(all(C ~= 0, 2), :);
    nf = size(C,1);
    nc = G.cells.num;
    D = sparse([(1:nf)'; (1:nf)'], C, ones(nf,1)*[-1 1], nf, nc);
    grad = @(x) D*x;
    div = @(x) -D'*x;
    avg = @(x) 0.5 * (x(C(:,1)) + x(C(:,2)));
    
    % Residual Calculation
    K=exp(-avg(T));
    res = div(K.*grad(T));
    
    % Imposing Boundary Condition
    % left cells: res = T - 1
    res(left_boundary) = T(left_boundary) - 1;
    % right cells: res = T - 0
    res(right_boundary) = T(right_boundary) - 0;
end