% mainGenSequences();
% load(fullfile('data', 'dummyDNA.mat'));
% pcPWMp = misc.preComputePWMp(X);
% mainPWM(pcPWMp, X, Y);
% mergedPeaksMin = mainGenSequences(1000, 600, 2, true);
% mergedPeaksMin = load(fullfile('data', 'dummyDNA.mat'));
% cd /cs/stud/boogalla/cbioDavid/projects/CompGenetics/BaumWelch/src
%
% mergedPeaksMin = load('data/peaks/roadmap/mergedPeaksMinimized.mat');
% mergedPeaksMin = mergedPeaksMin.mergedPeaksMin;
% mainRealData(mergedPeaksMin);


function mainRealData(mergedPeaksMin)
    dbstop if error
    close all;
    params.m = 1;
    params.order = 3;
    [params.k, params.n, params.J] = size(misc.PWMs());
    % params.tEpsilon = 1 ./ (mergedPeaksMin.lengths(1));
    params.tEpsilon = 0;
    params.batchSize = 2;
    testTrainRatio = 0.10;

    [test, train] = preprocess(mergedPeaksMin, testTrainRatio);
    % rocAucTest(params, train.pcPWMp, train.Y);

    learnedThetas = {};
    realM = max(train.Y, [] , 1);
    for i = 1:realM
        % train each base state
        X = train.X(train.Y(:, 1, 1)==i, :);
        pcPWMp = train.pcPWMp(train.Y(:, 1, 1)==i, :, :);
        params.m = 1;
        learnedThetas{i} = learnSingleMode(X, params, pcPWMp, 2);
        theta = learnedThetas{i};
    end
    learnedTheta = catThetas(params, learnedThetas);
    params.m = realM;
    accuricy = classify(learnedTheta, params, train.X, train.pcPWMp, train.Y)
    title('Train Dataset');
    accuricy = classify(learnedTheta, params, test.X, test.pcPWMp, test.Y)
    title('Test Dataset');


    % % merge thetas
    % testParams = params;
    % testParams.m = length(unique(test.Y));
    % testTheta = misc.genThetaJ(testParams);
    % for i=1:max(train.Y, [], 1);
    %     testTheta.E(i,:) = thetas(i).E(:);
    %     testTheta.G(i,:) = thetas(i).G(:);
    %     testTheta.F(i) = thetas(i).F;
    % end
    % % test
    % % N x m x L + J
    % figure
    % subplot(1, 2, 1);
    % scatter(1:length(mergedPeaksMin.originalTheta.E(:)), exp(mergedPeaksMin.originalTheta.E(:)))
    % hold on;
    % scatter(1:length(testTheta.E(:)), exp(testTheta.E(:)));
    % legend('Original Theta','Trained Theta')
    % subplot(1, 2, 2);
    % scatter(1:length(mergedPeaksMin.originalTheta.G(:)), exp(mergedPeaksMin.originalTheta.G(:)))
    % hold on;
    % scatter(1:length(testTheta.G(:)), exp(testTheta.G(:)));
    % legend('Original Theta','Trained Theta')

    % t1 = mergedPeaksMin.originalTheta;
    % t2 = testTheta;
    % alpha1 = EM.forwardAlgJ(train.X, t1, testParams, train.pcPWMp);
    % beta1 = EM.backwardAlgJ(train.X, t1, testParams, train.pcPWMp);
    % pX1 = EM.makePx(alpha1, beta1);
    % alpha2 = EM.forwardAlgJ(train.X, t2, testParams, train.pcPWMp);
    % beta2 = EM.backwardAlgJ(train.X, t2, testParams, train.pcPWMp);
    % pX2 = EM.makePx(alpha2, beta2);
    % % t2.E = t1.E;
    % figure
    % subplot(1,4,1);
    % hold on
    % plot(permute(alpha1(1,1,:), [3,2,1]))
    % plot(permute(alpha2(1,1,:), [3,2,1]))
    % legend('Original Theta','Trained Theta')
    % subplot(1,4,2);
    % hold on
    % plot(permute(beta1(1,1,:), [3,2,1]))
    % plot(permute(beta2(1,1,:), [3,2,1]))
    % legend('Original Theta','Trained Theta')
    % subplot(1,4,3);
    % hold on
    % scatter(1:length(pX1), pX1);
    % scatter(1:length(pX2), pX2);
    % legend('Original Theta','Trained Theta')
    % subplot(1,4,4);
    % plot(permute(alpha1(1,1,:)+beta1(1,1,:), [3,2,1]))
    % hold on
    % plot(permute(alpha2(1,1,:)+beta2(1,1,:), [3,2,1]))
    % legend('Original Theta','Trained Theta');
    % keyboard
    % [~, YsEst] = max(theta.gamma(:,:,1:end-params.J), [], 10);
    % YsEst = permute(YsEst, [1,3,2]);
    % calcError(Y(:)', YsEst(:)');
end

function rocAucTest(params, pcPWMp, Y)

    sortedPWMs = sort(pcPWMp, 3);
    PWMmax = sum(sortedPWMs(:, :, end-10: end), 3);
    pos = PWMmax(Y == 1, :);
    neg = PWMmax(Y ~= 1, :);
    aucRocs = zeros(params.k, 1);
    aucRocsSign = false(params.k, 1);
    for i = 1 : params.k
        [aucRocs(i),aucRocsSign(i)] = matUtils.getAucRoc(pos(:, i), neg(:, i), false, true);
        fprintf('%d - %.2f\n', i, aucRocs(i))
    end
    [PWM, lengths, names] = misc.PWMs();
    [b, i] = max(aucRocs, [], 1);

    fprintf('The best: %s - %.2f\n', names{i}, b)

    [bests, is] = sort(aucRocs, 1);
    sortedPWMs(:, aucRocsSign==1, :) = -sortedPWMs(:, aucRocsSign==1, :);
    PWMmaxOfBest = sum(sum(sortedPWMs(:, is(end-2:end), end-10 : end), 3), 2);
    pos = PWMmaxOfBest(Y == 1);
    neg = PWMmaxOfBest(Y ~= 1);
    fprintf('sum of best: %.2f\n', matUtils.getAucRoc(pos, neg, false, true))
end
% gamma - N x m x L
% psi - N x m x k x L
function Yest = genEstimation(params, theta, gamma, psi)
    [N, ~, L] = size(gamma);
    % N x L x m
    gammaPer = permute(gamma, [1, 3, 2]);
    [gammaMaxVals, YestBaseStates] = max(gammaPer, [], 3);

    % N x m x k x L -> N x L x m x k
    psiPer = cat(2, -inf(N, params.J, params.m, params.k), permute(psi, [1, 4, 2, 3]));
    skewedPsi = -inf(N, L, params.m, params.k);
    for l = 1:params.k
        for u = 1:theta.lengths(l)
            skewedPsi(:, :, :, l) = max(skewedPsi(:, :, :, l), psiPer(:, [1:L] + params.J - u, :, l));
        end
    end
    % N x L x m x k -> N x L
    [psiMaxVals, psiMaxInd] = max(skewedPsi(:,:,:), [], 3);
    subStates = floor((psiMaxInd - 1) / params.m) + 1;
    baseStates = mod(psiMaxInd - 1, params.m) + 1;
    subStateMask = psiMaxVals > gammaMaxVals;
    % N x L
    YestBaseStates(subStateMask) = baseStates(subStateMask);
    YestSubStates = zeros(N, L);
    YestSubStates(subStateMask) = subStates(subStateMask);
    Yest = cat(3, YestBaseStates, YestSubStates);

    % Yest1 = max(matUtils.logAdd(matUtils.logMatSum(skewedPsi, 4), gammaPer), [], 3);

end

function accuricy = classify(theta, params, X, pcPWMp, Y)
    [N, L] = size(X);
    fprintf('Calculating alpha...\n')
    % N x m x L
    alpha = EM.forwardAlgJ(X, theta, params, pcPWMp);
    fprintf('Calculating beta...\n')
    beta = EM.backwardAlgJ(X, theta, params, pcPWMp);
    % N x 1
    pX = EM.makePx(alpha, beta);
    fprintf('Calculating Gamma...\n')
    % N x m x L
    gamma = EM.makeGamma(params, alpha, beta, pX);
    % N x m x k x L
    psi = EM.makePsi(alpha, beta, X, params, theta, pcPWMp, pX);
    EM.drawStatus(theta, params, gamma)
    % g = permute(gamma, [1, 3, 2]);
    % figure;
    % plot(g(1, :, 1));
    % hold on;
    % plot(g(1, :, 2));
    % xlim([1,L])
    % % YPer = permute(Y, [3,2,1]);
    % % imagesc(YPer(:, :, 1)); colorbar;
    % legend('1', '2')
    % y=get(gca,'YLim');
    % for i = find(Y(1, :, 1) == 1)
    %     patch('XData', [i, i + 1, i + 1, i], 'YData',[y(1), y(1), y(2) y(2)],...
    %      'FaceColor', 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    % end


    psiSub = psi - repmat(pX, [1, params.m, params.k, L]);
    Yest = genEstimation(params, theta, gamma, psiSub);
    YestSingle = mode(Yest(:, :, 1), 2);
    Ymatch = Y == YestSingle;

    % figure;
    % subplot(1,3,1);
    % Ypresent = Y(:, :, 1);
    % Ypresent(~YmatchBoth) = params.m + 1;
    % imagesc(Ypresent);
    % colormap([1,0,0; 1,1,1; 1,0.5,0.5]);
    % caxis([1,3])
    % title('State Estimation')
    % drawnow;

    % subplot(1,3,2);
    % Ypresent = Y(:, :, 1);
    % Ypresent(~YmatchBase) = params.m + 1;
    % imagesc(Ypresent);
    % colormap([1,0,0; 1,1,1; 1,0.5,0.5]);
    % caxis([1,3])
    % drawnow;
    % title('Floor Estimation')

    % errorLengths = [];
    % isInError = false;
    % for i = 1:N
    %     for j = 1:L
    %         if YmatchBase(i, j) == true
    %             isInError = false;
    %         else
    %             if isInError
    %                 errorLengths(end) = errorLengths(end) + 1;
    %             else
    %                 errorLengths(end + 1) = 1;
    %             end
    %             isInError = true;
    %         end
    %     end
    %     isInError = false;
    % end
    % subplot(1,3,3);
    % histogram(errorLengths, 100, 'Normalization', 'probability');
    % title('Floor Errors Lengths Distribution')
    % mean(errorLengths)

    accuricy = mean(Ymatch(:), 1);
end

function [test, train] = preprocess(mergedPeaksMin, testTrainRatio)
    L = size(mergedPeaksMin.seqs, 2);
    types = [14, 15];
    overlaps = mergedPeaksMin.overlaps(:, :);
    % overlaps = mergedPeaksMin.overlaps(:, [2, 30]);
    mask = true(size(overlaps, 1), 1);
    mask = mask & mergedPeaksMin.lengths <= L;
    mask = mask & (sum(overlaps > 0, 2) == 1);
    mask = mask & (sum(overlaps(:, types) > 0, 2) == 1);
    % mask = mask & mergedPeaksMin.Y(:,1,1) == 1;
    mask = mask & mod(1:size(mask,1), 5).' == 0;
    overlaps = overlaps(mask, :);
    overlaps = overlaps(:, types);
    X = mergedPeaksMin.seqs(mask, :);
    [overlaps, seqInd] = sortrows(overlaps);
    X = X(seqInd, :);
    Y = (overlaps > 0) * (1:size(overlaps, 2))';

    % X = cat(2, X, fliplr(5-X));
    % N x k x L
    pcPWMp = misc.preComputePWMp(X);
    N = size(X, 1);
    trainMask = rand(N, 1) > testTrainRatio;
    train.X = X(trainMask, :);
    train.Y = Y(trainMask);
    train.pcPWMp = pcPWMp(trainMask, :, :);
    test.X = X(~trainMask, :);
    test.Y = Y(~trainMask);
    test.pcPWMp = pcPWMp(~trainMask, :, :);
end


% pcPWMp - N x k x L-1+J
% X - N x L

function [theta] = learnSingleMode(X, params, pcPWMp, maxIter)
    [theta, ~] = EM.EMJ(X, params, pcPWMp, maxIter);
end

function theta = meanMergeTheta(params, thetas)
    thetas = [thetas{:}];
    parts = length(thetas);
    theta = misc.genThetaUni(params);
    theta.G = log(mean(reshape(exp([thetas.G]), [params.m, params.k, parts]), 3));
    theta.F = log(mean(exp([thetas.F]), 2));
    theta.T = log(mean(reshape(exp([thetas.T]), [params.m, params.m, parts]), 3));
    theta.E = zeros([params.m, ones(1, params.order) * params.n]);
    for i = 1:parts
        theta.E = theta.E + exp(thetas(i).E);
    end
    theta.E = log(theta.E ./ parts);
end

function theta = catThetas(params, thetas)
    thetas = [thetas{:}];
    params.m = length(thetas);
    theta = misc.genThetaUni(params);
    theta.G = reshape([thetas.G], [params.k, params.m])';
    %theta.F = [thetas.F]';
    theta.T = log(eye(params.m) * (1 - (params.m * params.tEpsilon)) + params.tEpsilon);
    theta.E = zeros([params.m, ones(1, params.order) * params.n]);
    for i = 1:params.m
        theta.E(i, :) = thetas(i).E(:);
    end
end

% P - u x 1
% Q - u x 1
function ret = relativeEntropy(P, Q)
    ret = sum(P .* (log(P) - log(Q)), 1);
end

% % P - u x 1
% % Q - u x 1
% function ret = mse(P, Q)
%     ret = mean((P - Q) .^ 2, 1);
% end