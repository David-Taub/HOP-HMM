
% seqs - N x L
% indices - (L - order + 1) * N  x 1 (numbers from 1 to order) the indices in E
function indices = getIndices1D(seqs, order)
    [N, L] = size(seqs);
    matSize = 4 * ones(1, order);

    k = zeros(N, L - order + 1, order);
    for i = 1 : order
        k(:, :, i) = seqs(:, i : end - order + i);
    end
    % order x L - order + 1 x N
    k = permute(k, [3, 2, 1]);
    % order x (L - order + 1) * N -> (L - order + 1) * N  x 1
    indices = matUtils.matSub2ind(matSize, k(:, :));
end

