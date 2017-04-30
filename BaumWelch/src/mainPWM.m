% mainGenSequences();
% load(fullfile('data', 'dummyDNA.mat'));
% pcPWMp = BaumWelchPWM.preComputePWMp(X);
% mainPWM(pcPWMp, X, Y);

% pcPWMp - N x k x L-1
function mainPWM(pcPWMp, Xs, Ys)
    params.m = 5;
    params.order = 3;
    params.n = max(Xs(:));
    [params.N, params.L] = size(Xs);
    params.J = size(pcPWMp, 3) - params.L + 1;
    params.k = size(pcPWMp, 2);
    params.tEpsilon = 0.01;
    [theta, ~] = learn(Xs, params, pcPWMp);
    % N x m x L + J
    [~, YsEst] = max(theta.gamma(:,:,1:end-params.J), [], 10);
    YsEst = permute(YsEst, [1,3,2]);
    keyboard
    calcError(Ys(:)', YsEst(:)')
end

% L - sequence lengths
% means unique overlap between tissues, and if n is 1:3 then we take
% the sequences of the three most frequent class
function [posSeqs, negSeqs] = loadTommySeqs(L)
    fprintf('Loading sequences of Tommy...\n');
    negSeqsTrain = matUtils.readSeq(fullfile('data', 'NEnhancers.train.seq'), L);
    posSeqsTrain = matUtils.readSeq(fullfile('data', 'Enhancers.train.seq'), L);
    negSeqsTest = matUtils.readSeq(fullfile('data', 'NEnhancers.test.seq'), L);
    posSeqsTest = matUtils.readSeq(fullfile('data', 'Enhancers.test.seq'), L);

    posSeqs = [posSeqsTest; posSeqsTrain];
    negSeqs = [negSeqsTest; negSeqsTrain];

    parmas.N = min(size(negSeqs, 1), size(posSeqs, 1));
    posSeqs = posSeqs(1:parmas.N, :);
    negSeqs = negSeqs(1:parmas.N, :);
    % shuffle
    negSeqs = negSeqs(randperm(size(negSeqs, 1)), :);
    posSeqs = posSeqs(randperm(size(posSeqs, 1)), :);
end

% pcPWMp - N x k x L-1+J
% XTrain - N x L
function [theta, likelihood] = learn(Xs, params, pcPWMp)
    maxIter = 20;
    [theta, likelihood] = BaumWelchPWM.EMJ(Xs, params, pcPWMp, maxIter);
end

