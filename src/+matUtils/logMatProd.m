function ret = logMatProd(A, B)
    [n1, ~] = size(A);
    [~, n2] = size(B);
    % n1 x n3 x n2
    Arep = repmat(A, [1, 1, n2]);
    % n1 x n3 x n2
    Brep = repmat(permute(B, [3, 1, 2]), [n1, 1, 1]);
    ret = matUtils.logMatSum(Arep + Brep, 2);
    ret = permute(ret, [1, 3, 2]);
end