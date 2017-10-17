% alpha - N x m x L
% beta - N x m x L
% psi - N x m x k x L
% pX - N x 1
% gamma - N x m x L
function gamma = makeGamma(params, alpha, beta, pX)
    gamma = alpha + beta;
    gamma = gamma - repmat(pX, [1, params.m, size(alpha, 3)]);
end
