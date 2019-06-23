function XIndicesHotMap =  genXInidcesHotMap(params, dataset)
    [N, L] = size(dataset.X);
    % N x L - order + 1
    indices = reshape(matUtils.getIndices1D(dataset.X, params.order, params.n), [L - params.order + 1, N]).';
    % N x L - order + 1 x maxEIndex
    XIndicesHotMap = matUtils.mat23Dmat(indices, params.n ^ params.order);
    % N x L  x maxEIndex
    XIndicesHotMap = cat(2, false(N, params.order - 1, params.n ^ params.order), XIndicesHotMap);
end
