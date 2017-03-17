% res - N x k
function res = getPWMp(PWMstep, Xs1H, t)
    J = size(PWMstep, 2);
    % N x J x n
    lastJXs1H = Xs1H(:, t+1:t+J, :);
    % N x J x n x k, N x J x n
    res =  sumDim(bsxfun(@times, PWMsRep, lastJXs1H), [2, 3]);
end