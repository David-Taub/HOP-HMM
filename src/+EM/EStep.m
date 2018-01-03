% X - N x L
% alpha - N x m x L
% beta - N x m x L
% pX - N x 1
% xi - N x m x m x L
% gamma - N x m x L
% psi - N x m x k x L
function [alpha, beta, pX, xi, gamma, psi] = EStep(params, theta, X, pcPWMp)
    % fprintf('Running E step...\n');
    [N, L] = size(X);
    Eps = EM.getEp3d(theta, params, X, 1:L);
    alpha = EM.forwardAlg(X, theta, params, pcPWMp, Eps);
    % fprintf('Calculating beta...\n')
    beta = EM.backwardAlg(X, theta, params, pcPWMp, Eps);
    % N x 1
    pX = EM.makePx(alpha, beta);
    % fprintf('Calculating Xi...\n')
    % xi - N x m x m x L
    xi = EM.makeXi(theta, params, alpha, beta, X, pX);
    % fprintf('Calculating Gamma...\n')
    % gamma - N x m x L
    gamma = EM.makeGamma(params, alpha, beta, pX);
    % N x m x k x L
    psi = EM.makePsi(alpha, beta, X, params, theta, pcPWMp, pX);

end