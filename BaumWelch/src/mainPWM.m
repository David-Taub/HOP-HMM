function mainPWM()
    % get the n'st most frequent overlap
    L = 500;
    n = [1:1];
    [posSeqs, negSeqs] = loadTommySeqs(n, L);

    trainLabLength = ceil(size(posSeqs,1) / 2);

    XTrain = [posSeqs(1:trainLabLength, :); negSeqs(1:trainLabLength, :)];
    XTest = [posSeqs(trainLabLength + 1: end, :); negSeqs(trainLabLength + 1:end, :)];
    YTrain = cat(1, ones(trainLabLength, 1), ones(trainLabLength, 1) .* 2);
    YTest = cat(1, ones(size(posSeqs, 1) - trainLabLength,1), ones(size(posSeqs, 1) - trainLabLength,1) .* 2);
    learn(XTrain);


end

% L - sequence lengths
% n - number of classes to get for the positive sequences, where class
% means unique overlap between tissues, and if n is 1:3 then we take 
% the sequences of the three most frequent class
function [posSeqs, negSeqs] = loadTommySeqs(n, L)
    negSeqsTrain = matUtils.readSeq(fullfile('Data', 'NEnhancers.train.seq'), L);
    posSeqsTrain = matUtils.readSeq(fullfile('Data', 'Enhancers.train.seq'), L);
    negSeqsTest = matUtils.readSeq(fullfile('Data', 'NEnhancers.test.seq'), L);
    posSeqsTest = matUtils.readSeq(fullfile('Data', 'Enhancers.test.seq'), L);

    posSeqs = [posSeqsTest; posSeqsTrain];
    negSeqs = [negSeqsTest; negSeqsTrain];

    N = min(size(negSeqs, 1), size(posSeqs, 1));
    posSeqs = posSeqs(1:N, :);
    negSeqs = negSeqs(1:N, :);
    % shuffle
    negSeqs = negSeqs(randperm(size(negSeqs, 1)), :);
    posSeqs = posSeqs(randperm(size(posSeqs, 1)), :);
end

function [datasets] = learn(posSeqs)
    order = 3; m = 10; n = 4; maxIter = 300; tEpsilon = 0.01;
    [PWMs, lengths] = BaumWelchPWM.PWMs();
    [startT, T, Y, E, F, likelihood, gamma] = BaumWelchPWM.EMJ(posSeqs, m, maxIter, tEpsilon, order, PWMs, lengths)
end
