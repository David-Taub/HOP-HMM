
% seqs - N x L
% indices - 1 x n (numbers from 1 to order)
function indices = getIndeices1D(seqs, order)
    [N, L] = size(seqs);
    matSize = 4 * ones(1, order);

    k = zeros(N, L - order + 1, order);
    for i = 1 : order
        k(:, :, i) = seqs(:, i : end - order + i);
    end
    k = permute(k, [3, 1, 2]);
    indices = matSub2ind(matSize, k(:, :));
end

