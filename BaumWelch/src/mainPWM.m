% mainGenSequences();
% load(fullfile('data', 'dummyDNA.mat'));
% pcPWMp = BaumWelchPWM.preComputePWMp(X);
% mainPWM(pcPWMp, X, Y);
% cd /cs/stud/boogalla/cbioDavid/projects/CompGenetics/BaumWelch/src
% mergedPeaksMin = load('data/peaks/roadmap/mergedPeaksMinimized.mat');
% mainGenSequences(20, 50);
% mergedPeaksMin = load('data/dummyDNA.mat');
% mainPWM(mergedPeaksMin);


function mainPWM(mergedPeaksMin)
    close all;
    [test, train] = preprocess(mergedPeaksMin);
    trainParams.m = 1;
    trainParams.order = 3;
    trainParams.n = max(train.X(:));
    [trainParams.N, trainParams.L] = size(train.X);
    trainParams.J = size(BaumWelchPWM.PWMs(), 3);
    trainParams.k = size(train.pcPWMp, 2);
    trainParams.tEpsilon = 0.0001;

    for i = 1:max(train.Y)
        % train each base state
        X = train.X(train.Y==i, :);
        pcPWMp = train.pcPWMp(train.Y==i, :, :);
        trainParams.N = sum(train.Y==i, 1);
        [theta, ~] = learn(X, trainParams, pcPWMp, mergedPeaksMin.originalTheta, 3);
        thetas(i) = theta;
    end
    % merge thetas
    testParams = trainParams;
    testParams.m = length(unique(test.Y));
    testParams.N = size(test.X, 1);
    testTheta = BaumWelchPWM.genThetaJ(testParams);
    for i=unique(train.Y);
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
    subplot(1,2,2);
    scatter(1:length(testTheta.G(:)), testTheta.G(:));
    hold on;
    scatter(1:length(mergedPeaksMin.originalTheta.G(:)), mergedPeaksMin.originalTheta.G(:))

    t1 = mergedPeaksMin.originalTheta;
    t2 = testTheta;
    alpha1 = BaumWelchPWM.EM.forwardAlgJ(train.X, t1, trainParams, train.pcPWMp);
    beta1 = BaumWelchPWM.EM.backwardAlgJ(train.X, t1, trainParams, train.pcPWMp);
    pX1 = BaumWelchPWM.EM.makePx(alpha1, beta1);
    alpha2 = BaumWelchPWM.EM.forwardAlgJ(train.X, t2, trainParams, train.pcPWMp);
    beta2 = BaumWelchPWM.EM.backwardAlgJ(train.X, t2, trainParams, train.pcPWMp);
    pX2 = BaumWelchPWM.EM.makePx(alpha2, beta2);
    % t2.E = t1.E;
    % alpha3 = BaumWelchPWM.EM.forwardAlgJ(train.X, t2, trainParams, train.pcPWMp);
    % beta3 = BaumWelchPWM.EM.backwardAlgJ(train.X, t2, trainParams, train.pcPWMp);
    % pX3 = BaumWelchPWM.EM.makePx(alpha3, beta3);
    figure
    subplot(1,4,1);
    hold on
    plot(permute(alpha1(1,1,:), [3,2,1]))
    plot(permute(alpha2(1,1,:), [3,2,1]))
    % plot(permute(alpha3(1,1,:), [3,2,1]))
    legend('Original Theta','Trained Theta')
    subplot(1,4,2);
    hold on
    plot(permute(beta1(1,1,:), [3,2,1]))
    plot(permute(beta2(1,1,:), [3,2,1]))
    % plot(permute(beta3(1,1,:), [3,2,1]))
    legend('Original Theta','Trained Theta')
    % legend('1','2','3')
    subplot(1,4,3);
    hold on
    scatter(1:length(pX1), pX1);
    scatter(length(pX1)+1:length(pX1)*2, pX2);
    legend('Original Theta','Trained Theta')
    subplot(1,4,4);
    plot(permute(alpha1(1,1,:)+beta1(1,1,:), [3,2,1]))
    plot(permute(alpha2(1,1,:)+beta2(1,1,:), [3,2,1]))
    legend('Original Theta','Trained Theta')
    classify(testTheta, testParams, test.X, test.pcPWMp, test.Y);
    % [~, YsEst] = max(theta.gamma(:,:,1:end-params.J), [], 10);
    % YsEst = permute(YsEst, [1,3,2]);
    % calcError(Y(:)', YsEst(:)');
end

function classify(theta, params, X, pcPWMp, Y)
    fprintf('Calculating alpha...\n')
    alpha = BaumWelchPWM.EM.forwardAlgJ(X, theta, params, pcPWMp);
    fprintf('Calculating beta...\n')
    beta = BaumWelchPWM.EM.backwardAlgJ(X, theta, params, pcPWMp);
    % N x 1
    pX = BaumWelchPWM.EM.makePx(alpha, beta);
    fprintf('Calculating Gamma...\n')
    % gamma - N x m x L
    gamma = BaumWelchPWM.EM.makeGamma(params, alpha, beta, pX);
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

