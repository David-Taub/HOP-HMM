% mainGenSequences();
% load(fullfile('data', 'dummyDNA.mat'));
% pcPWMp = BaumWelchPWM.preComputePWMp(X);
% mainPWM(pcPWMp, X, Y);
% cd /cs/stud/boogalla/cbioDavid/projects/CompGenetics/BaumWelch/src
% mergedPeaksMin = load('data/peaks/roadmap/mergedPeaksMinimized.mat');
% mergedPeaksMin = mainGenSequences(70, 90, 1);
% mainPWM(mergedPeaksMin);


function mainPWM(mergedPeaksMin)
    close all;
    [test, train] = preprocess(mergedPeaksMin);
    trainParams.m = 1;
    trainParams.order = 3;
    trainParams.n = max(train.X(:));
    trainParams.J = size(BaumWelchPWM.PWMs(), 3);
    trainParams.k = size(train.pcPWMp, 2);
    trainParams.tEpsilon = 0.0001;

    for i = 1:max(train.Y)
        % train each base state
        X = train.X(train.Y==i, :);
        pcPWMp = train.pcPWMp(train.Y==i, :, :);
        trainParams.m = 1;
        trainTheta = BaumWelchPWM.genThetaUni(trainParams);
        % trainTheta.G = mergedPeaksMin.originalTheta.G;
        [theta, ~] = learn(X, trainParams, pcPWMp, trainTheta, 2);
        thetas(i) = theta;
    end
    % merge thetas
    testParams = trainParams;
    testParams.m = length(unique(test.Y));
    testTheta = BaumWelchPWM.genThetaJ(testParams);
    for i=1:max(train.Y, [], 1);
        testTheta.E(i,:) = thetas(i).E(:);
        testTheta.G(i,:) = thetas(i).G(:);
        testTheta.F(i) = thetas(i).F;
    end
    % test
    % N x m x L + J
    figure
    subplot(1, 2, 1);
    scatter(1:length(testTheta.E(:)), testTheta.E(:));
    hold on;
    scatter(1:length(mergedPeaksMin.originalTheta.E(:)), mergedPeaksMin.originalTheta.E(:))
    subplot(1, 2, 2);
    scatter(1:length(testTheta.G(:)), testTheta.G(:));
    hold on;
    scatter(1:length(mergedPeaksMin.originalTheta.G(:)), mergedPeaksMin.originalTheta.G(:))

    t1 = mergedPeaksMin.originalTheta;
    t2 = testTheta;
    alpha1 = BaumWelchPWM.EM.forwardAlgJ(train.X, t1, testParams, train.pcPWMp);
    beta1 = BaumWelchPWM.EM.backwardAlgJ(train.X, t1, testParams, train.pcPWMp);
    pX1 = BaumWelchPWM.EM.makePx(alpha1, beta1);
    alpha2 = BaumWelchPWM.EM.forwardAlgJ(train.X, t2, testParams, train.pcPWMp);
    beta2 = BaumWelchPWM.EM.backwardAlgJ(train.X, t2, testParams, train.pcPWMp);
    pX2 = BaumWelchPWM.EM.makePx(alpha2, beta2);
    % t2.E = t1.E;
    figure
    subplot(1,4,1);
    hold on
    plot(permute(alpha1(1,1,:), [3,2,1]))
    plot(permute(alpha2(1,1,:), [3,2,1]))
    legend('Original Theta','Trained Theta')
    subplot(1,4,2);
    hold on
    plot(permute(beta1(1,1,:), [3,2,1]))
    plot(permute(beta2(1,1,:), [3,2,1]))
    legend('Original Theta','Trained Theta')
    subplot(1,4,3);
    hold on
    scatter(1:length(pX1), pX1);
    scatter(1:length(pX2), pX2);
    legend('Original Theta','Trained Theta')
    subplot(1,4,4);
    plot(permute(alpha1(1,1,:)+beta1(1,1,:), [3,2,1]))
    hold on
    plot(permute(alpha2(1,1,:)+beta2(1,1,:), [3,2,1]))
    legend('Original Theta','Trained Theta');
    classify(testTheta, testParams, test.X, test.pcPWMp, test.Y)
    classify(testTheta, testParams, train.X, train.pcPWMp, train.Y)
    keyboard
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
    figure
    hold on
    [~, YEst] = max(gamma(:, :, 1), [], 2);
    accuracy = sum(Y == YEst) / length(Y);
end

function [test, train] = preprocess(mergedPeaksMin)
    L = size(mergedPeaksMin.seqs, 2);
    overlaps = mergedPeaksMin.overlaps(:, :);
    % overlaps = mergedPeaksMin.overlaps(:, [1,2,3,4]);
    mask = mergedPeaksMin.lengths >= L;
    mask = mask & (sum(overlaps > 0, 2) == 1);
    % mask = mask & mod(1:size(mask,1), 15).' == 0;
    overlaps = overlaps(mask, :);
    X = mergedPeaksMin.seqs(mask, :);
    Y = (overlaps > 0) * (1:size(overlaps, 2))';
    [overlaps, seqInd] = sortrows(overlaps);
    X = X(seqInd, :);

    % X = cat(2, X, fliplr(5-X));
    fprintf('Calculating PWMs LogLikelihood\n')
    size(X)
    % N x k x L
    pcPWMp = BaumWelchPWM.preComputePWMp(X);
    N = size(X, 1);
    trainMask = rand(N, 1) > 0.1;
    train.X = X(trainMask, :);
    train.Y = Y(trainMask);
    train.pcPWMp = pcPWMp(trainMask, :, :);
    test.X = X(~trainMask, :);
    test.Y = Y(~trainMask);
    test.pcPWMp = pcPWMp(~trainMask, :, :);
end


% pcPWMp - N x k x L-1+J
% XTrain - N x L
function [theta, likelihood] = learn(X, params, pcPWMp, initTheta, maxIter)
    [theta, likelihood] = BaumWelchPWM.EM.EMJ(X, params, pcPWMp, initTheta, maxIter);
end

