
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
    indicesHotMap = cat(2, false(params.N, params.order-1+params.J, params.n ^ params.order), indicesHotMap);
    theta = BaumWelchPWM.genThetaJ(params);
    % drawStatus(theta, params, 0,0,0);
t    iterLike = [];
    for it = 1:maxIter;
        tic
        % alphaBase - N x m x L
        % alphaSub - N x m x L+J x k
        % N x m x L
        alpha = BaumWelchPWM.forwardAlgJ(Xs, theta, params, pcPWMp);
        beta = BaumWelchPWM.backwardAlgJ(Xs, theta, params, pcPWMp);
        pX = makePx(theta, params, alpha, beta);
        xi = makeXi(theta, params, alpha, beta, Xs, pX);
        gamma = makeGamma(theta, params, alpha, beta, pX);
        psi = makePsi(theta, params, alpha, beta, pcPWMpRep, Xs);

        assert(not(any(isnan(theta.T(:)))))
        assert(not(any(isnan(theta.E(:)))))
        assert(not(any(isnan(theta.G(:)))))
        assert(not(any(isnan(theta.F(:)))))
        assert(not(any(isnan(alpha(:)))))
        assert(not(any(isnan(beta(:)))))
        theta.E = updateE(gamma, theta, params, indicesHotMap);
        theta.T = updateT(xi, gamma, params);
        theta.startT = updateStartT(gamma);
        theta.G = updateG(psi);
        % [theta, gamma] = updateTheta(theta, params, Xs, indicesHotMap, pcPWMp, alphaBase, alphaSub, beta);
        % iterLike(end + 1) = sum(log(scale(:)));
        iterLike(end+1) = sum(pX, 1)
        % DRAW
        drawStatus(theta, params, alpha, beta, gamma);


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

% psi - N x m x k x L
% pcPWMp - N x k x L
% alpha - N x m x L
% beta - N x m x L
function psi = makePsi(theta, params, alpha, beta, pcPWMpRep, Xs)
    kronMN = kron(1:params.m, ones(1, params.N));
    matSize = [params.m , params.n * ones(1, params.order)];
    % Ret - N x m x L
    Eps = BaumWelchPWM.getEp3d(theta, params, Xs, 1:params.L, kronMN, matSize);

    psi = repmat(permute(alpha, [1, 2, 4, 3]), [1, 1, params.k, 1]);
    psi  = psi + repmat(theta.F', [params.N, 1, params.k, params.L]);
    psi  = psi + repmat(permute(theta.G, [3, 1, 2]), [params.N, 1, 1, params.L]);
    psi  = psi + repmat(permute(Eps, [1, 2, 4, 3]), [1, 1, params.k, 1]);
    for l = 1:params.k
        psi(:, :, l, :) = psi(:, :, l, :) + cat(3, beta(:, :, params.lengths(l):end), -inf(params.N, params.m, params.lengths(l)));
    end
end


% alpha - N x m x L
% beta - N x m x L
% psi - N x m x k x L
% pX - N x 1
% gamma - N x m x L
function gamma = makeGamma(theta, params, alpha, beta, pX)
    gamma = alpha + beta;
    gamma = logMatSum(gamma, 2);
    gamma = gamma - repmat(pX, [1, params.m, params.m, params.L]);
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
    xi  = xi + repmat(compF, [params.N, 1, params.m, L]);
    xi  = xi + repmat(permute(theta.T, [3, 1, 2]), [params.N, 1, 1, params.L]);
    xi  = xi + repmat(permute(theta.G, [3, 1, 2]), [params.N, 1, 1, params.L]);
    xi  = xi + repmat(permute(Eps, [1, 4, 2, 3]), [1, params.m, 1, 1]);
    betaMoved = cat(3, beta(:, :, 2:end), -inf(params.N, params.m, 1))
    xi  = xi + repmat(permute(betaMoved, [1, 4, 2, 3]), [1, params.m, 1, 1]);
    xi  = xi - repmat(pX, [1, params.m, params.m, params.L]);
end

% pX - N x 1
function pX = makePx(theta, params, alpha, beta)
    % N x m x L
    pX = alpha + beta;
    pX = logMatSum(pX, 2);
    pX = pX(:, 1, 1)
end

function drawStatus(theta, params, alpha, beta, gamma)
    [~, YsEst] = max(gamma(:,:,1:end), [], 2);
    subplot(2,3,1);imagesc(permute(YsEst, [1,3,2])); colorbar;
    title('gamma')
    [~, YsEst] = max(alpha(:,:,1:end), [], 2);
    subplot(2,3,2);imagesc(permute(YsEst, [1,3,2])); colorbar;
    title('alpha')
    [~, YsEst] = max(beta(:,:,1:end), [], 2);
    subplot(2,3,3);imagesc(permute(YsEst, [1,3,2])); colorbar;
    title('beta')
    subplot(2,3,4);imagesc(theta.F); colorbar;
    title('F')
    subplot(2,3,5);imagesc(theta.T); colorbar;
    title('T')
    subplot(2,3,6);imagesc(theta.G); colorbar;
    title('M')
    drawnow
    % keyboard
end

function [theta, gamma] = updateTheta(theta, params, Xs, indicesHotMap, pcPWMp, alphaBase, alphaSub, beta)
    fprintf('Updating theta\n');
    % N x m x L + J
    % gamma_t(i) = P(y_t = i|x_1:L)
    % gamma = alphaBase .* beta ./ (repmat(sum(alpha.Base .* beta, 2), [1, params.m, 1]) + eps);
    gamma = alphaBase .* beta;
    [theta.T, theta.G, theta.F] = updateTGF(theta, params, Xs, pcPWMp, alphaBase, alphaSub, beta);
    theta.E = updateE(gamma, params, indicesHotMap);
    theta.startT = updateStartT(gamma);
    theta.T = updateT(gamma, xi);
    % T bound trick
    theta.T = Tbound(theta, params);
end

% indicesHotMap - N x L x maxEIndex
% gamma - N x m x L
% E - m x 4 x 4 x 4 x ... x 4 (order times)
function newE = updateE(gamma, params, indicesHotMap)
    %  m x N x L + J
    fprintf('Update E\n');
    perGamma = permute(gamma, [2,1,3]);
    newE = -inf([params.m, params.n * ones(1, params.order)]);
    for i = 1:(params.n ^ params.order)
        % m x N x L -> m x 1
        newE(:, i) = sum(exp(perGamma(:, indicesHotMap(:,:,i))), 2);
    end
    newE = log(bsxfun(@times, newE, 1 ./ sum(newE, params.order+1)));
end

% gamma - N x m x L
function newStartT = updateStartT(gamma)
    newStartT = sum(exp(gamma(:,:,1)), 1);
    % probability distribution normalization
    newStartT = log(newStartT / sum(newStartT, 2)).';
end

% xi - N x m x m x L
% gamma - N x m x L
% newT - m x m
function newT = updateT(xi, gamma, params)
    newT = permute(sum(sum(exp(xi), 1), 4), [2,3,1]);
    newT = newT ./ repmat(sum(sum(exp(gamma), 1), 3), [params.m,1]).';
    newT = log(bsxfun(@times, newT, 1 ./ sum(newT, 2)));
end

% psi - N x m x k x L
% gamma - N x m x L
% newG - m x k
function newG = updateG(psi)
    newG = permute(sum(sum(exp(psi), 1), 4), [2,3,1]);
    newG = log(bsxfun(@times, newG, 1 ./ sum(newG, 2)));
end

% newT - m x m
function newT = Tbound(theta, params)
    x = 0.5 * params.tEpsilon * (params.m-1);
    for i = 1 : params.m
        theta.T(i, :) = (theta.T(i, :) * x) / ((sum(theta.T(i, :),2)-theta.T(i, i)) * (1 + x));
        theta.T(i, i) = 1 / (x + 1);
        % for j = 1 : params.m
        %     if i ~= j && theta.T(i,j) > params.tEpsilon
        %         theta.T(i, i) = theta.T(i, i) + (theta.T(i, j) - params.tEpsilon);
        %         theta.T(i, j) = params.tEpsilon;
        %     end
        % end
    end
    newT = theta.T;
end
