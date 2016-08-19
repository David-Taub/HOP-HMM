function main()
    close all;
    L = 500;
    posSeqsTrain = readSeq('Enhancers.train.seq', L);
    negSeqsTrain = readSeq('NEnhancers.train.seq', L);
    posSeqsTest = readSeq('Enhancers.test.seq', L);
    negSeqsTest = readSeq('NEnhancers.test.seq', L);

    Xtrain = cat(1, posSeqsTrain, negSeqsTrain);
    Xtest = cat(1, posSeqsTest, negSeqsTest);
    YTrain = cat(1, ones(size(posSeqsTrain, 1),1), zeros(size(negSeqsTrain, 1),1));
    YTest = cat(1, ones(size(posSeqsTest, 1),1), zeros(size(negSeqsTest, 1),1));

    m = 2;
    parfor order = 6:6
        
        [posE, negE] = trainMarkov(posSeqsTrain, negSeqsTrain, order);
        
        thresholds = 0.98 : 0.00005 : 1.01;
        [trainErr, threshold] = classify(posE, negE, posSeqsTrain, negSeqsTrain, thresholds);
        [testErr, threshold] = classify(posE, negE, posSeqsTest, negSeqsTest, threshold);
        [order, threshold, trainErr, testErr]

        pos2neg = 1 / 300;
        neg2pos = 1 / 50;
        
        [startT, T, E] = createHmmParams(posE, negE, neg2pos, pos2neg);

        % N x 1
        posPostirior = getPostirior(posSeqsTrain, startT, T, E);
        negPostirior = getPostirior(negSeqsTrain, startT, T, E);
        % N x 1
        posTops = getTopPart(posPostirior);
        negTops = getTopPart(negPostirior);

        minTops = min(min(posTops), min(negTops));
        maxTops = max(max(posTops), max(negTops));
        success = [];
        thresholds = minTops : 0.01 : maxTops;
        [trainErr, threshold] = findThreshold(posTops, negTops, thresholds);
        % N x 1
        posPostirior = getPostirior(posSeqsTest, startT, T, E);
        negPostirior = getPostirior(negSeqsTest, startT, T, E);
        % N x 1
        posTops = getTopPart(posPostirior);
        negTops = getTopPart(negPostirior);
        testErr = getLose(posTops, negTops, threshold);
        [order, trainErr, testErr]

        % figure 
        % hold on
        % plot(posTops)
        % plot(negTops)
        % hold off
        % legend('pos', 'neg')
        % figure
        % hold on;
        % plot(mean(posPostirior, 1));
        % plot(mean(negPostirior, 1));
        % ylim([0,1]);
        % legend('positive postirior', 'negative postirior');
        % title('postirior Probability of Being Enhancer');
        % hold off;
    end
    save('data.mat')
end


function out = getTopPart(M)
    L = size(M, 2);
    marginsRatio = 0.10;
    marginsRatio2 = 0.01;
    topPartRatio = 0.15;
    M = M(:, ceil(L * marginsRatio) : end - ceil(L * marginsRatio));
    M = sort(M, 2, 'descend');
    M = M(:, ceil(L * marginsRatio2) : end - ceil(L * marginsRatio2));
    out = mean(M(: ,1:ceil(L * topPartRatio)), 2);
end

% high - N1 x 1
% low - N2 x 1
% thresholds - 1 x R
function [err, threshold] = findThreshold(high, low, thresholds)
    N = size(high, 1) + size(low, 1);
    % N1 x R
    tp = bsxfun(@lt, repmat(high, [1, length(thresholds)]), thresholds);
    % N2 x R
    tn = bsxfun(@ge, repmat(low, [1, length(thresholds)]), thresholds);
    % 1 x R
    errs = (sum(tp, 1) + sum(tn, 1)) ./ N;
    [err, i] = min(errs, [], 2);
    threshold = thresholds(i);
end

% high - N1 x 1
% low - N2 x 1
% thresholds - 1 x R
function err = getLose(high, low, threshold)
    N = size(high, 1) + size(low, 1);
    err = (sum(high < threshold, 1) + sum(low >= threshold, 1)) ./ N;
end
% seqs - S x L
% N - number of sequences to calculate the postirior with
% out - N x L
function out = getPostirior(seqs, startT, T, E)
    m = 2;
    [N, L] = size(seqs);
    postirior = zeros(m, L, N);
    [alpha, scale] = forwardAlg(seqs, startT, T, E);
    beta = backwardAlg(seqs, startT, T, E, scale);
    % S x m x L
    postirior = alpha .* beta;
    postirior = bsxfun(@times, postirior, 1 ./ sum(postirior, 2));
    % return postirior of the positive state, 
    out(:,:) = postirior(:, 1, :);
end

function [startT, T, E] = createHmmParams(posE, negE, neg2pos, pos2neg)
    E = cat(1, shiftdim(posE, -1), shiftdim(negE, -1));
    T = [1 - pos2neg, pos2neg; neg2pos, 1 - neg2pos];
    startT = [0.5; 0.5];
end

function [posE, negE] = trainMarkov(posSeqs, negSeqs, order)
    posE = getEFromSeqs(posSeqs, order);
    negE = getEFromSeqs(negSeqs, order);
end

function [err, threshold] = classify(posE, negE, posSeqs, negSeqs, thresholds)
    likePosIsPos = getLogLikes(posE, posSeqs); %high
    likePosIsNeg = getLogLikes(negE, posSeqs); %low
    likeNegIsPos = getLogLikes(posE, negSeqs); %low
    likeNegIsNeg = getLogLikes(negE, negSeqs); %high

    ratioPos = likePosIsPos ./ likePosIsNeg; %high / low = high
    ratioNeg = likeNegIsPos ./ likeNegIsNeg; %low
    if length(thresholds) > 1
        [err, threshold] = findThreshold(ratioPos, ratioNeg, thresholds);
    else
        threshold = thresholds;
        err = getLose(ratioPos, ratioNeg, threshold);
    end
        
end

function E = getEFromSeqs(seqs, order)
    ambient = 10 ^ -6;
    N = length(seqs);
    matSize = [4 * ones(1, order), 1];
    E = zeros(matSize);
    for i = 1 : N
        indices = getIndeices1D(seqs(i, :), order);
        h = histc(indices, 1 : 4 ^ order);
        Ecur = reshape(h, [matSize, 1]);
        E = E + Ecur;
    end
    E = E + ambient;
    E = bsxfun(@times, E, 1 ./ sum(E, order));
end

function logLikes = getLogLikes(E, seqs)
    N = size(seqs, 1);
    logLikes = zeros(N,1);

    order = matDim(E);
    indices = getIndeices1D(seqs, order);

    STOPPED HERE
    STOPPED HERE
    STOPPED HERE
    STOPPED HERE
    STOPPED HERE
    STOPPED HERE
    indices = reshape(indices, [L, N]);
    logLike = sum(sum(log(E(indices)), 2), 1);
    STOPPED HERE
    STOPPED HERE
    STOPPED HERE
    STOPPED HERE
    STOPPED HERE
    for i=1:N
        logLikes(i) = getLogLikeFromSeq(seqs(i, :), E);
    end
end

% reads .seq format file (Tommy's format) and returns the sequence as numbers
function out = readSeq(filePath, L)
    % fprintf('Reading %s...', filePath);
    [~, seq] = textread(filePath, '%s%s%*[^\n]','delimiter','\t','bufsize',20000);
    out = zeros(length(seq), L); 
    for i=1:length(seq)
        seq{i}=upper(seq{i}); 
        if length(seq{i}) < L
            continue;
        end
        seq{i}=seq{i}(ceil(length(seq{i})/2) + [-floor(L/2) + 1: floor(L/2)]);
        out(i, :) = nt2int(seq{i});
    end
    out( ~any(out,2), : ) = [];  %remove zero rows
end


% seqs - N x L
% indices - 1 x n (numbers from 1 to order)
function indices = getIndeices1D(seqs, order)
    [N, L] = size(seqs);
    matSize = 4 * ones(1, order);

    k = zeros(N, L - order + 1, order);
    for i = 1 : order
        k(:, :, i) = seqs(:, i : end - order + i);
    end
    k = permute(k, [3, 2, 1]);

    indices = matSub2ind(matSize, k(:, :));
end

% seq - 1 x n
% E - 4 x ... x 4 (order times)
% logLike - number
function logLike = getLogLikeFromSeq(seqs, E)
    order = matDim(E);
    indices = getIndeices1D(seqs, order);
    logLike = sum(sum(log(E(indices)), 2), 1);
end
