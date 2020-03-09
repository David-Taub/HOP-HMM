% maxIter - maximal EM iterations allowed without convergence before quitting
% patience - maximal EM iterations allowed without likelihood improvement before quitting
% repeat - the amount of repeats with different theta initializations should be done before
%           returning the best performing theta
% In dataset:
%   X - N x L emission variables
%   pcPWMp - N x k x L
% In params:
%   m - amount of possible states (y)
%   n - amount of possible emissions (x)
%   order - the HMM order of the E matrix
%
function [bestTheta, bestLikelihood, bestThetas] = EM(dataset, params, maxIter, patience, repeat)

    [N, L] = size(dataset.X);
    fprintf('Starting EM algorithm on %d x %d\n', N, L);
    bestLikelihood = -Inf;
    % figure
    for rep = 1:repeat
        % X = X(randperm(N), :);
        fprintf('Repeat %d / %d\n', rep, repeat);
        initTheta = misc.genTheta(params, false, false);
        [iterationLikelihood, theta, thetas] = EMRun(dataset, params, initTheta, maxIter, params.doGTBound, patience);
        if bestLikelihood < iterationLikelihood
            bestLikelihood = iterationLikelihood;
            bestTheta = theta;
            bestThetas = thetas;
        end
    end
end



% iterates the EM algorithm, returns the likelihood of the best iteration, and theta parameters at that iteration
function [iterationLikelihood, theta, thetas] = EMRun(dataset, params, initTheta, maxIter, doGTBound, patience)
    LIKELIHOOD_THRESHOLD = 10 ^ -6;
    thetas(1) = initTheta;
    iterLikes = -inf(maxIter, 1);

    for iterationIndex = 1:maxIter
        fprintf('EM iteration %d / %d\n', iterationIndex, maxIter);
        tic;
        [thetas(iterationIndex + 1), iterationLikelihood] = EM.EMIteration(params, dataset, thetas(iterationIndex), doGTBound);
        motifsPer = sum(exp(thetas(iterationIndex + 1).G(:)), 1) .* 100;
        timeLapse = toc();
        fprintf('iterationIndex %d: log-like: %.2f Time: %.2fs, motifs: ~%.2f%%\n', iterationIndex, iterationLikelihood, timeLapse, motifsPer);
        if iterationIndex > patience && iterationIndex > 1 && ...
            abs((iterationLikelihood - iterLikes(iterationIndex - 1)) / iterationLikelihood) < LIKELIHOOD_THRESHOLD
            fprintf('Converged\n');
            break
        end

        if iterationIndex > patience && iterationLikelihood < iterLikes(iterationIndex - patience)
            fprintf('Patience reached, Converged\n');
            keyboard
            break
        end
        iterLikes(iterationIndex) = iterationLikelihood;
    end % end EM iteration loop
    [iterationLikelihood, ind] = max(iterLikes);
    theta = thetas(ind);
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