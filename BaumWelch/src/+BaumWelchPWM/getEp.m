function Ep = getEp(E, Xs, t, m, kronMN, matSize, N, order)
    if t < order
        E = sumDim(E, [2 : 1 + order - t]);
        subscripts = [kronMN; repmat(Xs(:, 1 : t).', [1, m])];
        indices = matSub2ind(matSize(1 : t + 1), subscripts);
    else
        subscripts = [kronMN; repmat(Xs(:, t-order+1 : t).', [1, m])];
        indices = matSub2ind(matSize, subscripts);
    end
    % N x m
    Ep = reshape(E(indices).', [N, m]);
end