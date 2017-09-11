% mainGenSequences();
% load(fullfile('data', 'dummyDNA.mat'));
% pcPWMp = BaumWelchPWM.preComputePWMp(X);
% mainPWM(pcPWMp, X, Y);
% cd /cs/stud/boogalla/cbioDavid/projects/CompGenetics/BaumWelch/src
% mergedPeaksMin = load('data/peaks/roadmap/mergedPeaksMinimized.mat');
% mergedPeaksMin = mainGenSequences(2000, 1000, 1);
% mainPWM(mergedPeaksMin);


function mainPWM(mergedPeaksMin)
    close all;
    trainParams.m = 1;
    trainParams.order = 3;
    [trainParams.k, trainParams.n, trainParams.J] = size(BaumWelchPWM.PWMs());
    trainParams.tEpsilon = 0.0001;
    testTrainRatio = 0.10;
    [test, train] = preprocess(mergedPeaksMin, testTrainRatio);
    learnedThetas = {};
    for i = 1:max(train.Y(:))
        % train each base state
        X = train.X(train.Y(:, 1)==i, :);
        pcPWMp = train.pcPWMp(train.Y(:, 1)==i, :, :);
        learnedThetas{i} = learnSingleMode(X, trainParams, pcPWMp, 3, 10, mergedPeaksMin);
    end
    learnedTheta = catThetas(trainParams.params, learnedThetas);
    classify(learnedTheta, testParams, test.X, test.pcPWMp, test.Y)
    classify(learnedTheta, testParams, train.X, train.pcPWMp, train.Y)



    % % merge thetas
    % testParams = trainParams;
    % testParams.m = length(unique(test.Y));
    % testTheta = BaumWelchPWM.genThetaJ(testParams);
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
    % alpha1 = BaumWelchPWM.EM.forwardAlgJ(train.X, t1, testParams, train.pcPWMp);
    % beta1 = BaumWelchPWM.EM.backwardAlgJ(train.X, t1, testParams, train.pcPWMp);
    % pX1 = BaumWelchPWM.EM.makePx(alpha1, beta1);
    % alpha2 = BaumWelchPWM.EM.forwardAlgJ(train.X, t2, testParams, train.pcPWMp);
    % beta2 = BaumWelchPWM.EM.backwardAlgJ(train.X, t2, testParams, train.pcPWMp);
    % pX2 = BaumWelchPWM.EM.makePx(alpha2, beta2);
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

function accuracy = classify(theta, params, X, pcPWMp, Y)
    fprintf('Calculating alpha...\n')
    alpha = BaumWelchPWM.EM.forwardAlgJ(X, theta, params, pcPWMp);
    fprintf('Calculating beta...\n')
    beta = BaumWelchPWM.EM.backwardAlgJ(X, theta, params, pcPWMp);
    % N x 1
    pX = BaumWelchPWM.EM.makePx(alpha, beta);
    fprintf('Calculating Gamma...\n')
    % gamma - N x m x L
    gamma = BaumWelchPWM.EM.makeGamma(params, alpha, beta, pX);
    % figure
    % hold on
    [~, YEst] = max(gamma(:, :, 1), [], 2);
    accuracy = sum(Y == YEst) / length(Y);
end

function [test, train] = preprocess(mergedPeaksMin, testTrainRatio)
    L = size(mergedPeaksMin.seqs, 2);
    overlaps = mergedPeaksMin.overlaps(:, :);
    % overlaps = mergedPeaksMin.overlaps(:, [1,2,3,4]);
    mask = mergedPeaksMin.lengths >= L;
    mask = mask & (sum(overlaps > 0, 2) == 1);
    % mask = mask & mod(1:size(mask,1), 15).' == 0;
    overlaps = overlaps(mask, :);
    X = mergedPeaksMin.seqs(mask, :);
    Y = mergedPeaksMin.Y(mask, :);
    % Y = (overlaps > 0) * (1:size(overlaps, 2))';
    [overlaps, seqInd] = sortrows(overlaps);
    X = X(seqInd, :);
    Y = Y(seqInd, :);

    % X = cat(2, X, fliplr(5-X));
    % N x k x L
    pcPWMp = BaumWelchPWM.preComputePWMp(X);
    N = size(X, 1);
    trainMask = rand(N, 1) > testTrainRatio;
    train.X = X(trainMask, :);
    train.Y = Y(trainMask, :);
    train.pcPWMp = pcPWMp(trainMask, :, :);
    test.X = X(~trainMask, :);
    test.Y = Y(~trainMask, :);
    test.pcPWMp = pcPWMp(~trainMask, :, :);
end


% pcPWMp - N x k x L-1+J
% X - N x L
function [theta] = learnSingleMode(X, params, pcPWMp, maxIter, parts, mergedPeaksMin)
    thetas = {};
    N = size(X,1);
    dists = [];
    for i = 1:parts
        partMask = mod(1:N, parts) == i-1;
        % [test, train] = preprocess(mergedPeaksMin, testTrainRatio);
        initTheta = BaumWelchPWM.genThetaUni(params);
        initTheta.ot = mergedPeaksMin.originalTheta;
        % theta.Y = train.Y;
        [thetas{i}, ~] = BaumWelchPWM.EM.EMJ(X(partMask, :), params, pcPWMp(partMask, :, :), initTheta, maxIter);

        theta = meanMergeTheta(params, thetas);
        dists(i) = relativeEntropy(exp(mergedPeaksMin.originalTheta.G'), exp(theta.G'));
        dists(i) = dists(i) + relativeEntropy(exp(mergedPeaksMin.originalTheta.E(:)), exp(theta.E(:)));
        figure
        subplot(1,3,1);
        plot(dists);
        title('KLdiv over iterations')
        subplot(1,3,2);
        title('G')
        plot(exp(mergedPeaksMin.originalTheta.G));ylim([0,1]);
        hold on;
        plot(exp(theta.G));ylim([0,1]);
        legend('Original Theta','Trained Theta')
        subplot(1,3,3);
        title('E')
        plot(exp(mergedPeaksMin.originalTheta.E(:)));ylim([0,1]);
        hold on;
        plot(exp(theta.E(:)));ylim([0,1]);
        legend('Original Theta','Trained Theta')
        drawnow
    end
    theta = meanMergeTheta(params, thetas);
end

function theta = meanMergeTheta(params, thetas)
    thetas = [thetas{:}];
    parts = length(thetas);
    theta = BaumWelchPWM.genThetaUni(params);
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
    theta = BaumWelchPWM.genThetaUni(params);
    theta.G = reshape([thetas.G], [params.k, params.m])';
    theta.F = [thetas.F]';
    theta.T = reshape([thetas.T], [params.m, params.m])';
    theta.E = zeros([params.m, ones(1,order) * params.n]);
    for i = 1:params.m
        theta.E(1, :) = thetas(i).E(:);
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