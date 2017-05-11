% res - N x k
function res = getPWMp(params, PWMsRep, Xs1H, t, mask)
    % N x J x n
    lastJXs1H = Xs1H(:,t:t+params.J-1,:);
    % N x J x n x k, N x J x n -> N x k
    res = matUtils.sumDim(bsxfun(@times, PWMsRep, lastJXs1H), 3);
    res = matUtils.mulDim(res + mask, 2);
end