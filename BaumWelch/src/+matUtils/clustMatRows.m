function [M, leafOrd] = clustMatRows(M)
    tree = linkage(M);
    leafOrd = optimalleaforder(tree, pdist(M));
    M = M(leafOrd, :);
end
