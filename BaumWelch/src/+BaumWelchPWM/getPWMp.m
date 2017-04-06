% res - N x k
function res = getPWMp(PWMsRep, Xs1H, t)
    J = size(PWMsRep, 2);
    % N x J x n
    lastJXs1H = Xs1H(:, t:t+J-1, :);
    % N x J x n x k, N x J x n -> N x k
    res =  matUtils.sumDim(bsxfun(@times, PWMsRep, lastJXs1H), [2, 3]);
    
end