% X - N x L emission variables
% m - amount of possible states (y)
% n - amount of possible emissions (x)
% maxIter - maximal iterations allowed
% tEpsilon - probability to switch states will not exceed this number (if tEpsilon = 0.01,
%            and m = 3 then the probability to stay in a state will not be less than 0.98)
% pcPWMp - N x k x L
% order - the HMM order of the E matrix
% initial estimation parameters
function [bestTheta, bestLikelihood] = EM(dataset, params, maxIter, patience, repeat)

    [N, L] = size(dataset.X);
    fprintf('Starting EM algorithm on %d x %d\n', N, L);
    bestLikelihood = -Inf;
    % figure
    for rep = 1:repeat
        % X = X(randperm(N), :);
        fprintf('Repeat %d / %d\n', rep, repeat);
        initTheta = misc.genTheta(params, true);
        [iterLike, theta] = EMRun(dataset, params, initTheta, maxIter, false, patience);
        if bestLikelihood < iterLike
            bestLikelihood = iterLike;
            bestTheta = theta;
        end
    end
    if params.doGTBound
        [bestLikelihood, bestTheta] = EMRun(dataset, params, initTheta, maxIter, true, patience);
    end
end



% iterates the EM algorithm, returns the likelihood of the best iteration, and theta parameters at that iteration
function [iterLike, theta] = EMRun(dataset, params, initTheta, maxIter, doGTBound, patience)
    LIKELIHOOD_THRESHOLD = 10 ^ -6;
    theta(1) = initTheta;
    iterLikes = -inf(maxIter, 1);

    for it = 1:maxIter
        fprintf('EM iteration %d / %d\n', it, maxIter);
        tic;
        [theta(it + 1), iterLike] = EM.EMIteration(params, dataset, theta(it), doGTBound);
        motifsPer = sum(exp(theta(it + 1).G(:)), 1).*100;
        timeLapse = toc();
        fprintf('It %d: log-like: %.2f Time: %.2fs, motifs: ~%.2f%%\n', it, iterLike, timeLapse, motifsPer);
        if it > 1 && abs((iterLike - iterLikes(it - 1)) / iterLike) < LIKELIHOOD_THRESHOLD
            fprintf('Converged\n');
            break
        end

        if it > patience && iterLike < iterLikes(it - patience)
            fprintf('Patience reached, Converged\n');
            break
        end
        iterLikes(it) = iterLike;
    end % end EM iteration loop
    [iterLike, ind] = max(iterLikes);
    theta = theta(ind);
 end


% psi - N x m x k x L
% pX - N x 1
% gamma - N x m x L
function drawStatus(theta, params, alpha, beta, gamma, pX, xi, psi)
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    YsEst = mean(gamma, 3);
    YsEst = YsEst(:, 1);
    subplot(2,2,1);
    scatter(1:size(pX, 1), YsEst); colorbar;
    title('Px');
    % [~, YsEst] = max(alpha, [], 2);
    % subplot(2,3,2);imagesc(permute(YsEst, [1,3,2])); colorbar;
    % % title('alpha')
    % [~, YsEst] = max(beta(:,:,1:end), [], 2);
    % subplot(2,3,3);imagesc(permute(YsEst, [1,3,2])); colorbar;
    % title('beta')
    subplot(2,2,2);plot(exp(theta.E(:)));
    plot(exp(theta.E(:)));
    ylim([0,1]);
    title('E');

    subplot(2,2,3);
    plot(permute(gamma(1,1,:), [3,2,1]));
    hold on;
    % N x m x k x L
    plot(permute(matUtils.logMatSum(psi(1,1,:,:), 3), [4,3,2,1]));
    title('Postirior');
    subplot(2,2,4);plot(exp(theta.G(:)));
    ylim([0,1]);
    title('G');
    drawnow;
end