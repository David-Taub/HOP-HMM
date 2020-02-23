% Ep - N x m
function Ep = getEp(theta, params, X, t, kronMN, matSize)
    if t < params.order
        E = matUtils.sumDim(exp(theta.E), 2 : 1 + params.order - t);
        E = log(bsxfun(@times, E, 1 ./ sum(E, length(size(E)))));
        subscripts = [kronMN; repmat(X(:, 1 : t).', [1, params.m])];
        indices = matUtils.matSub2ind(matSize(1 : t + 1), subscripts);
    else
        E = theta.E;
        subscripts = [kronMN; repmat(X(:, t-params.order+1 : t).', [1, params.m])];
        indices = matUtils.matSub2ind(matSize, subscripts);
    end
    % N x m
    Ep = reshape(E(indices).', [size(X, 1), params.m]);
    assert(not(any(isnan(Ep(:)))));
end