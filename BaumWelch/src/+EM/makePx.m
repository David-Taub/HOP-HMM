
% pX - N x 1
function pX = makePx(alpha, beta)
    % N x m x L
    pX = alpha + beta;
    pX = matUtils.logMatSum(pX, 2);
    pX = pX(:, 1, 1);
end
