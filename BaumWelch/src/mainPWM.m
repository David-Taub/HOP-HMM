% mainGenSequences();
% load(fullfile('data', 'dummyDNA.mat'));
% pcPWMp = BaumWelchPWM.preComputePWMp(X);
% mainPWM(pcPWMp, X);

% pcPWMp - N x k x L-1
function mainPWM(pcPWMp, Xs)
    m = 5;order = 3;
    [~, lengths] = BaumWelchPWM.PWMs();
    learn(Xs, m, order, pcPWMp, lengths);

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

    N = min(size(negSeqs, 1), size(posSeqs, 1));
    posSeqs = posSeqs(1:N, :);
    negSeqs = negSeqs(1:N, :);
    % shuffle
    negSeqs = negSeqs(randperm(size(negSeqs, 1)), :);
    posSeqs = posSeqs(randperm(size(posSeqs, 1)), :);
end

% pcPWMp - N x k x L-1+J
% XTrain - N x L
function [startT, T, M, E, F, likelihood, gamma] = learn(Xs, m, order, pcPWMp, lengths)
    maxIter = 20; tEpsilon = 0.01;
    [startT, T, M, E, F, likelihood, gamma] = BaumWelchPWM.EMJ(Xs, m, order, pcPWMp, lengths, maxIter, tEpsilon);
end

