
function [bestTheta, bestLikelihood] = EMJ(Xs, params, pcPWMp, maxIter)
    % PWMs - k x n x J emission matrix of m Jaspar PWMs with length J
    %        true length of i'th PWM< J and given in lengths(i) if a PWM is
    %        shorter than j, it is aligned to the end of the 3rd dimension.
    % Xs - N x L emission variables
    % m - amount of possible states (y)
    % n - amount of possible emissions (x)
    % maxIter - maximal iterations allowed
    % tEpsilon - probability to switch states will not exceed this number (if tEpsilon = 0.01,
    %            and m = 3 then the probability to stay in a state will not be less than 0.98)
    % order - the HMM order of the E matrix
    % initial estimation parameters

    LIKELIHOOD_THRESHOLD = 10 ^ -4;
    bestLikelihood = -Inf;
    repeat = 1;
    pcPWMpRep = repmat(permute(pcPWMp, [1, 4, 2, 3]), [1, params.m, 1, 1]);

    % N x L -order + 1
    indices = reshape(matUtils.getIndices1D(Xs, params.order), [params.L-params.order+1, params.N]).';
    % N x L -order + 1 x maxEIndex
    indicesHotMap = matUtils.mat23Dmat(indices, params.n ^ params.order);
    % N x L  x maxEIndex
    indicesHotMap = cat(2, false(params.N, params.order-1, params.n ^ params.order), indicesHotMap);
    theta = BaumWelchPWM.genThetaJ(params);
    % drawStatus(theta, params, 0,0,0);
    iterLike = [];
    for it = 1:maxIter;
        tic
        % alphaBase - N x m x L
        % alphaSub - N x m x L+J x k
        % N x m x L
        fprintf('Calculating alpha...\n')
        alpha = BaumWelchPWM.forwardAlgJ(Xs, theta, params, pcPWMp);
        fprintf('Calculating beta...\n')
        beta = BaumWelchPWM.backwardAlgJ(Xs, theta, params, pcPWMp);
        % N x 1
        pX = makePx(alpha, beta);
        fprintf('Calculating Xi...\n')
        xi = makeXi(theta, params, alpha, beta, Xs, pX);
        fprintf('Calculating Gamma...\n')
        gamma = makeGamma(params, alpha, beta, pX);
        fprintf('Calculating Psi...\n')
        % psi = makePsi(theta, params, alpha, beta, Xs);

        theta.E = updateE(gamma, params, indicesHotMap);
        theta.T = updateT(xi, gamma, params);

        theta.startT = updateStartT(gamma);
        % theta.G = updateG(psi);
        theta.G = updateG(alpha, beta, Xs, params, theta);
        % [theta, gamma] = updateTheta(theta, params, Xs, indicesHotMap, pcPWMp, alphaBase, alphaSub, beta);
        % iterLike(end + 1) = sum(log(scale(:)));
        iterLike(end+1) = matUtils.logMatSum(pX, 1);
        % DRAW
        drawStatus(theta, params, alpha, beta, gamma);
        assert(not(any(isnan(theta.T(:)))))
        assert(not(any(isnan(theta.E(:)))))
        assert(not(any(isnan(theta.G(:)))))
        assert(not(any(isnan(theta.F(:)))))
        assert(not(any(isnan(alpha(:)))))
        assert(not(any(isnan(beta(:)))))


        fprintf('Likelihood in iteration %d is %.2f (%.2f seconds)\n', length(iterLike), iterLike(end), toc());
        if length(iterLike) > 1 && abs((iterLike(end) - iterLike(end-1)) / iterLike(end)) < LIKELIHOOD_THRESHOLD
            fprintf('Converged\n');
            break
        end

        if bestLikelihood < iterLike(end)
            bestLikelihood = iterLike(end);
            bestTheta = theta;
            bestTheta.gamma = gamma;
        end
    end % end of iterations loop
end

% alpha - N x m x L
% beta - N x m x L
% psi - N x m x k x L
function psi = makePsi(theta, params, alpha, beta, Xs)
    kronMN = kron(1:params.m, ones(1, params.N));
    matSize = [params.m , params.n * ones(1, params.order)];
    % Ret - N x m x L
    Eps = BaumWelchPWM.getEp3d(theta, params, Xs, 1:params.L, kronMN, matSize);

    psi = repmat(permute(alpha, [1, 2, 4, 3]), [1, 1, params.k, 1]);
    psi  = psi + repmat(theta.F', [params.N, 1, params.k, params.L]);
    psi  = psi + repmat(permute(theta.G, [3, 1, 2]), [params.N, 1, 1, params.L]);
    psi  = psi + repmat(permute(Eps, [1, 2, 4, 3]), [1, 1, params.k, 1]);
    for l = 1:params.k
        psi(:, :, l, :) = psi(:, :, l, :) + permute(cat(3, beta(:, :, theta.lengths(l)+1:end), -inf(params.N, params.m, theta.lengths(l))), [1,2,4,3]);
    end
end


% alpha - N x m x L
% beta - N x m x L
% psi - N x m x k x L
% pX - N x 1
% gamma - N x m x L
function gamma = makeGamma(params, alpha, beta, pX)
    gamma = alpha + beta;
    gamma = gamma - repmat(pX, [1, params.m, params.L]);
end


% alpha - N x m x L
% beta - N x m x L
% pX - N x 1
% xi - N x m x m x L
function xi = makeXi(theta, params, alpha, beta, Xs, pX)
    kronMN = kron(1:params.m, ones(1, params.N));
    matSize = [params.m , params.n * ones(1, params.order)];
    % Eps - N x m x L
    Eps = BaumWelchPWM.getEp3d(theta, params, Xs, 1:params.L, kronMN, matSize);
    compF = log(1-exp(theta.F));
    % xi - N x m x m x L
    xi = repmat(permute(alpha, [1, 2, 4, 3]), [1, 1, params.m, 1]);

    xi  = xi + repmat(compF', [params.N, 1, params.m, params.L]);
    xi  = xi + repmat(permute(theta.T, [3, 1, 2]), [params.N, 1, 1, params.L]);
    xi  = xi + repmat(permute(Eps, [1, 4, 2, 3]), [1, params.m, 1, 1]);
    betaMoved = cat(3, beta(:, :, 2:end), -inf(params.N, params.m, 1));
    xi  = xi + repmat(permute(betaMoved, [1, 4, 2, 3]), [1, params.m, 1, 1]);
    xi  = xi - repmat(pX, [1, params.m, params.m, params.L]);

end

% pX - N x 1
function pX = makePx(alpha, beta)
    % N x m x L
    pX = alpha + beta;
    pX = matUtils.logMatSum(pX, 2);
    pX = pX(:, 1, 1);
end

% psi - N x m x k x L
% pX - N x 1
% gamma - N x m x L
function drawStatus(theta, params, alpha, beta, gamma)
    YsEst = mean(gamma, 3);
    YsEst = YsEst(:, 1);
    subplot(2,2,1); scatter(1:params.N, YsEst); colorbar;
    title('gamma')
    % [~, YsEst] = max(alpha, [], 2);
    % subplot(2,3,2);imagesc(permute(YsEst, [1,3,2])); colorbar;
    % % title('alpha')
    % [~, YsEst] = max(beta(:,:,1:end), [], 2);
    % subplot(2,3,3);imagesc(permute(YsEst, [1,3,2])); colorbar;
    % title('beta')
    subplot(2,2,2);imagesc(theta.F); colorbar;
    title('F')
    subplot(2,2,3);imagesc(theta.T); colorbar;
    title('T')
    subplot(2,2,4);imagesc(theta.G); colorbar;
    title('M')
    drawnow
end

% function [theta, gamma] = updateTheta(theta, params, Xs, indicesHotMap, pcPWMp, alphaBase, alphaSub, beta)
%     fprintf('Updating theta\n');
%     % N x m x L + J
%     % gamma_t(i) = P(y_t = i|x_1:L)
%     % gamma = alphaBase .* beta ./ (repmat(sum(alpha.Base .* beta, 2), [1, params.m, 1]) + eps);
%     gamma = alphaBase .* beta;
%     [theta.T, theta.G, theta.F] = updateTGF(theta, params, Xs, pcPWMp, alphaBase, alphaSub, beta);
%     theta.E = updateE(gamma, params, indicesHotMap);
%     theta.startT = updateStartT(gamma);
%     theta.T = updateT(gamma, xi);
%     % T bound trick
%     theta.T = Tbound(theta, params);
% end

% indicesHotMap - N x L + J x maxEIndex
% gamma - N x m x L
% E - m x 4 x 4 x 4 x ... x 4 (order times)
function newE = updateE(gamma, params, indicesHotMap)
    %  m x N x L + J
    fprintf('Update E\n');
    % m x N x L
    perGamma = permute(gamma, [2,1,3]);
    newE = -inf([params.m, params.n * ones(1, params.order)]);
    for i = 1:(params.n ^ params.order)
        % m x N x L -> m x 1
        newE(:, i) = matUtils.logMatSum(perGamma(:, indicesHotMap(:,:,i)), 2);
    end
    newE = exp(newE);
    newE = log(bsxfun(@times, newE, 1 ./ sum(newE, params.order+1)));
end

% gamma - N x m x L
function newStartT = updateStartT(gamma)
    newStartT = exp(matUtils.logMatSum(gamma(:,:,1), 1));
    % probability distribution normalization
    newStartT = log(newStartT / sum(newStartT, 2)).';
end

% xi - N x m x m x L
% gamma - N x m x L
% newT - m x m
function newT = updateT(xi, gamma, params)
    newT = permute(matUtils.logMatSum(matUtils.logMatSum(xi, 1), 4), [2,3,1]);
    newT = newT - repmat(matUtils.logMatSum(matUtils.logMatSum(gamma, 1), 3), [params.m,1]).';
    newT = exp(newT);
    newT = bsxfun(@times, newT, 1 ./ sum(newT, 2));
    newT = log(Tbound(params, newT));
end

function newG = updateG(alpha, beta, Xs, params, theta)
    kronMN = kron(1:params.m, ones(1, params.N));
    matSize = [params.m , params.n * ones(1, params.order)];
    beta = cat(3, beta, -inf(params.N, params.m, params.J + 1));
    newG = -inf(1, params.m, params.k);
    for t = 1:params.L
        % N x m x k
        Ep = BaumWelchPWM.getEp(theta, params, Xs, t, kronMN, matSize);
        newPsi = repmat(alpha(:,:,t), [1, 1, params.k]);
        newPsi = newPsi + repmat(theta.F', [params.N, 1, params.k]);
        newPsi = newPsi + repmat(permute(theta.G, [3, 1, 2]), [params.N, 1, 1]);
        newPsi = newPsi + repmat(Ep, [1, 1, params.k]);
        for l = 1:params.k
            newPsi(:, :, l) = newPsi(:, :, l) + beta(:, :, t+theta.lengths(l)+1);
        end
        newG = matUtils.logAdd(newG, matUtils.logMatSum(newPsi, 1));
    end
    newG = exp(permute(newG, [2,3,1]));
    newG = log(bsxfun(@times, newG, 1 ./ sum(newG, 2)));
end
% psi - N x m x k x L
% gamma - N x m x L
% newG - m x k
% function newG = updateG(psi)
%     newG = permute(sum(sum(exp(psi), 1), 4), [2,3,1]);
%     newG = log(bsxfun(@times, newG, 1 ./ sum(newG, 2)));
% end

% newT - m x m
function newT = Tbound(params, T)
    for i = 1 : params.m
        if T(i, i) < 1-params.tEpsilon;
            T(i, :) = T(i, :) * (params.tEpsilon / (1-T(i, i)));
            T(i, i) = 1-params.tEpsilon;
        end
    end
    newT = T;
end
