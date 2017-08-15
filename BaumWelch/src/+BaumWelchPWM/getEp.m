function Ep = getEp(theta, params, Xs, t, kronMN, matSize)
    if t < params.order
        E = matUtils.logMatSumDims(theta.E, 2 : 1 + params.order - t);
        subscripts = [kronMN; repmat(Xs(:, 1 : t).', [1, params.m])];
        indices = matUtils.matSub2ind(matSize(1 : t + 1), subscripts);
    else
        E = theta.E;
        subscripts = [kronMN; repmat(Xs(:, t-params.order+1 : t).', [1, params.m])];
        indices = matUtils.matSub2ind(matSize, subscripts);
    end
    % N x m
    Ep = reshape(E(indices).', [params.N, params.m]);
end