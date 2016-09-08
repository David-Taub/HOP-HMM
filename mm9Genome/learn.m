function [MMmean, MMvar, amounts] = learn(posSeqs, negSeqs, overlaps)
    % get the n'st most frequent overlap
    % load('/cs/stud/boogalla/  projects/CompGenetics/BaumWelch/peaks.mat');
    % posSeqs = seqs;
    % negSeqs = readSeq('NEnhancers.seq', L);
    M = size(overlaps, 2); %23
    MMmean = zeros(1, M+1);
    MMvar = zeros(1, M+1);
    amounts = zeros(1, M+1);
    for overlapClass = 0:M
        [MMmean(overlapClass+1), MMvar(overlapClass+1), ~, ~, ~, amounts(overlapClass+1)] = sampleAndLearnMulti(posSeqs, negSeqs, overlaps, overlapClass);
    end
    % plot(mean(testErrMMs(1:k,:), 2));
    % plotEs(Es, M);
end

function [MMmean, MMvar, HMMmean, HMMvar, Emean, amount] = sampleAndLearnMulti(posSeqs, negSeqs, overlaps, overlapClass)

    order = 6;
    repeats = 5;

    testErrMMs  = zeros(1, repeats);
    testErrHMMs = zeros(1, repeats);
    Es = zeros(repeats, 4 ^ order);
    for i = 1:repeats
        [XTrain, XTest, YTrain, YTest] = loadSeqs(posSeqs, negSeqs, overlaps, overlapClass);
        amount = sum([YTrain == 1;YTest == 1], 1);
        [testErrMMs(i), testErrHMMs(i), E] = learnData(XTrain, YTrain, XTest, YTest, order);
        Es(i, :) = E(1, :);
    end
    MMmean = mean(testErrMMs, 2);
    MMvar = var(testErrMMs, 0, 2);
    HMMmean = mean(testErrHMMs, 2);
    HMMvar = var(testErrHMMs, 0, 2);
    Emean = mean(Es, 1);
end

function plotEs(Es, M)
    tissues = {'BAT', 'BMDM', 'BoneMarrow',...
               'CH12', 'Cerebellum', 'Cortex',...
               'E14', 'Heart-E14.5', 'Heart',...
               'Kidney', 'Limb-E14.5', 'Liver-E14.5',...
               'Liver', 'MEF', 'MEL', 'OlfactBulb',...
               'Placenta', 'SmIntestine', 'Spleen',...
               'Testis', 'Thymus', 'WholeBrain-E14.5', 'mESC'};

    Es = mean(Es, 2);
    EsRep = repmat(Es, [1, M, 1]);
    D = median((EsRep - permute(EsRep, [2,1,3])) .^ 2, 3);
    figure;
    imagesc(D);colorbar;
    title('Emission Diff Between Tissues');
    % set(gca,'YLim',[0 M],'YTick',1:12,'YTickLabel',months)
    ax = gca;
    ax.XLim = [0 M];
    ax.YLim = [0 M];
    ax.XTick = 1:M;
    ax.YTick = 1:M;
    ax.YTickLabel = tissues;
    ax.XTickLabel = tissues;
    ax.XTickLabelRotation=45;
    % set(gca,'YLim',[0 M],'YTick',1:M, 'XLim',[0 M],'XTick',1:M,...
    %         'YTickLabel', tissues, 'XTickLabel', tissues);
end

% L - sequence lengths
% n - number of classes to get for the positive sequences, where class
% means unique overlap between tissues, and if n is 1:3 then we take 
% the sequences of the three most frequent class
function [XTrain, XTest, YTrain, YTest] = loadSeqs(posSeqs, negSeqs, overlaps, overlapClass)
    if overlapClass > 0
        posSeqs = posSeqs(overlaps(:, overlapClass) == 1, :);
    end

    % posSeqs = regularSeqs(posSeqs, L);


    % TODO: my BG seqs are bad. why?
    % load('/cs/stud/boogalla/projects/CompGenetics/BaumWelch/bg.mat');
    % negSeqs = seqs(:, ceil(size(seqs, 2)/2) + [-floor(L/2) + 1: floor(L/2)]);
    % negSeqs = readSeq('NEnhancers.seq', L);
    % assumes more negative than positives
    N = min(size(negSeqs, 1), size(posSeqs, 1));
    trainTestRate = 0.9;
    trainLabLength = ceil(N * trainTestRate);

    % shuffle
    posSeqs = posSeqs(randperm(size(posSeqs, 1), N), :);
    negSeqs = negSeqs(randperm(size(negSeqs, 1), N), :);

    % get train dataset
    XTrain = [posSeqs(1:trainLabLength, :); negSeqs(1:trainLabLength, :)];
    YTrain = [ones(trainLabLength, 1); ones(trainLabLength, 1) .* 2];
    

    % get test dataset
    XTest  = [posSeqs(trainLabLength + 1: N, :); negSeqs(trainLabLength + 1:N, :)];
    YTest  = [ones(N - trainLabLength,1); ones(N - trainLabLength,1) .* 2];
end

% L - sequence lengths
% n - number of classes to get for the positive sequences, where class
% means unique overlap between tissues, and if n is 1:3 then we take 
% the sequences of the three most frequent class
function [XTrain, XTest, YTrain, YTest] = loadSeqs2(posSeqs, negSeqs, overlaps, overlapClass)
    if overlapClass > 0
        posSeqs = posSeqs(overlaps(:, overlapClass) == 1, :);
    end

    % posSeqs = regularSeqs(posSeqs, L);


    % TODO: my BG seqs are bad. why?
    % load('/cs/stud/boogalla/projects/CompGenetics/BaumWelch/bg.mat');
    % negSeqs = seqs(:, ceil(size(seqs, 2)/2) + [-floor(L/2) + 1: floor(L/2)]);
    % negSeqs = readSeq('NEnhancers.seq', L);
    
    trainTestRate = 0.9;
    trainLabLength = ceil(size(posSeqs,1) * trainTestRate);
        
    N = min(size(negSeqs, 1), size(posSeqs, 1));
    posSeqs = posSeqs(1:N, :);
    negSeqs = negSeqs(1:N, :);
    % shuffle
    negSeqs = negSeqs(randperm(size(negSeqs, 1)), :);
    posSeqs = posSeqs(randperm(size(posSeqs, 1)), :);


    XTrain = [posSeqs(1:trainLabLength, :); negSeqs(1:trainLabLength, :)];
    XTest = [posSeqs(trainLabLength + 1: end, :); negSeqs(trainLabLength + 1:end, :)];
    YTrain = cat(1, ones(trainLabLength, 1), ones(trainLabLength, 1) .* 2);
    YTest = cat(1, ones(size(posSeqs, 1) - trainLabLength,1), ones(size(posSeqs, 1) - trainLabLength,1) .* 2);
end


% L - sequence lengths
% n - number of classes to get for the positive sequences, where class
% means unique overlap between tissues, and if n is 1:3 then we take 
% the sequences of the three most frequent class
% function [posSeqs, negSeqs] = loadTommySeqs(L)
%     posSeqsTrain = readSeq('Enhancers.train.seq', L);
%     posSeqsTest = readSeq('Enhancers.test.seq', L);
%     negSeqs = readSeq('NEnhancers.seq', L);
%     % negSeqsTrain = readSeq('NEnhancers.train.seq', L);
%     % negSeqsTest = readSeq('NEnhancers.test.seq', L);

%     posSeqs = [posSeqsTest; posSeqsTrain];
%     % negSeqs = [negSeqsTest; negSeqsTrain];
%     [posSeqs, negSeqs] = suffleAndTrim(posSeqs, negSeqs);

% end

function [posSeqs, negSeqs] = suffleAndTrim(posSeqs, negSeqs)
end 

function [testErrMM, testErrHMM, E] = learnData(XTrain, YTrain, XTest, YTest, order)
    % anaFreq(XTrain(YTrain == 1, :), XTrain(YTrain == 2, :), order);
    E = trainMarkov(XTrain, YTrain, order);
    thresholds = 0.0 : 0.005 : 2;
    [~, threshold] = classify(E, XTrain, YTrain, thresholds);
    [testErrMM, ~] = classify(E, XTest, YTest, threshold);
    % [order, threshold, testErrMM, testErr]
    testErrHMM = 0;

    % pos2neg = 1 / 250; % this values minimizes training error
    % neg2pos = 1 / 50;
    
    % [startT, T] = createHmmParams(neg2pos, pos2neg);

    % % N x 1
    % posPostirior = getPostirior(XTrain(YTrain == 1, :), startT, T, E);
    % negPostirior = getPostirior(XTrain(YTrain == 2, :), startT, T, E);
    % % N x 1
    % posTops = getTopPart(posPostirior);
    % negTops = getTopPart(negPostirior);

    % minTops = min(min(posTops), min(negTops));
    % maxTops = max(max(posTops), max(negTops));
    % success = [];
    % thresholds = minTops : 0.01 : maxTops;
    % [trainErr, threshold] = findThreshold(posTops, negTops, thresholds);
    % % N x 1
    % posPostirior = getPostirior(XTest(YTest == 1, :), startT, T, E);
    % negPostirior = getPostirior(XTest(YTest == 2, :), startT, T, E);
    % % N x 1
    % posTops = getTopPart(posPostirior);
    % negTops = getTopPart(negPostirior);
    % testErrHMM = getLose(posTops, negTops, threshold);
    % [order, trainErr, testErrHMM]
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
    % save('data.mat')
end


% function out = getTopPart(M)
%     L = size(M, 2);
%     marginsRatio = 0.15;
%     marginsRatio2 = 0.10;
%     topPartRatio = 0.3;
%     M = M(:, ceil(L * marginsRatio) : end - ceil(L * marginsRatio));
%     M = sort(M, 2, 'descend');
%     M = M(:, ceil(L * marginsRatio2) : end - ceil(L * marginsRatio2));
%     out = mean(M(: ,1:ceil(L * topPartRatio)), 2);
% end

% high - N1 x 1
% low - N2 x 1
% thresholds - 1 x R
function [err, threshold] = findThreshold(high, low, thresholds)
    % figure

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
% function out = getPostirior(seqs, startT, T, E)
%     m = 2;
%     [N, L] = size(seqs);
%     postirior = zeros(m, L, N);
%     [alpha, scale] = forwardAlg(seqs, startT, T, E);
%     beta = backwardAlg(seqs, startT, T, E, scale);
%     % S x m x L
%     postirior = alpha .* beta;
%     postirior = bsxfun(@times, postirior, 1 ./ sum(postirior, 2));
%     % return postirior of the positive state, 
%     out(:,:) = postirior(:, 1, :);
% end

% function [startT, T] = createHmmParams(neg2pos, pos2neg)
%     T = [1 - pos2neg, pos2neg; neg2pos, 1 - neg2pos];
%     startT = [0.5; 0.5];
% end

function E = trainMarkov(X, Y, order)
    E = [];
    % for i = [unique(Y)]
    for i = 1:max(Y, [], 1)
        Ei = getEFromSeqs(X(Y == i, :), order);
        Ei = shiftdim(Ei, -1);
        E = cat(1, E, Ei);
    end
end

% 2 label classify using the log liklihood ratio
function [err, threshold] = classify(E, X, Y, thresholds)
    s = size(E);
    posE = reshape(E(1,:), s(2:end));
    negE = reshape(E(2,:), s(2:end));
    likePosIsPos = getLogLikes(posE, X(Y == 1, :)); %high
    likePosIsNeg = getLogLikes(negE, X(Y == 1, :)); %low
    ratioPos = likePosIsPos ./ likePosIsNeg; %high / low = high
    
    likeNegIsPos = getLogLikes(posE, X(Y == 2, :)); %low
    likeNegIsNeg = getLogLikes(negE, X(Y == 2, :)); %high
    ratioNeg = likeNegIsPos ./ likeNegIsNeg; %low


    if length(thresholds) > 1
        [err, threshold] = findThreshold( ratioNeg.', ratioPos.', thresholds);
    else
        threshold = thresholds;
        err = getLose(ratioNeg.', ratioPos.', threshold);
    end
    
end

function E = getEFromSeqs(seqs, order)
    % ambient is a trick to avoid zero division for absent motifs
    ambient = 10 ^ -6;
    matSize = [4 * ones(1, order), 1];
    indices = getIndeices1D(seqs, order);
    h = histc(indices, 1 : 4 ^ order);
    E = reshape(h, [matSize, 1]) + ambient;
    E = bsxfun(@times, E, 1 ./ sum(E, order));
end

function anaFreq(seqsPos, seqsNeg, order)
    matSize = [4 * ones(1, order), 1];
    indicesP = getIndeices1D(seqsPos, order);
    indicesN = getIndeices1D(seqsNeg, order);
    hP = histc(indicesP, 1 : 4 ^ order);
    hN = histc(indicesN, 1 : 4 ^ order);
    h = hP - hN;
    N = 10;
    [sortedX,sortingIndices] = sort(h,'descend');
    % maxValues = sortedX(1:N);
    maxValueIndices = sortingIndices(1:N);
    maxValues = h(maxValueIndices);
    if order == 7
    [a(:, 1), a(:, 2), a(:, 3), a(:, 4), a(:, 5), a(:, 6), a(:, 7)] = ...
        ind2sub(matSize, maxValueIndices);
    elseif order == 6
    [a(:, 1), a(:, 2), a(:, 3), a(:, 4), a(:, 5), a(:, 6)] = ...
        ind2sub(matSize, maxValueIndices);
    elseif order == 5
    [a(:, 1), a(:, 2), a(:, 3), a(:, 4), a(:, 5)] = ...
        ind2sub(matSize, maxValueIndices);
    elseif order == 4
    [a(:, 1), a(:, 2), a(:, 3), a(:, 4)] = ...
        ind2sub(matSize, maxValueIndices);
    end
    for i = 1:N
        fprintf('%s (%d)\n', int2nt(a(i, :)), maxValues(i));
    end
    fprintf('\n');
    plot(sortedX)
    hold on
    drawnow

end


function logLikes = getLogLikes(E, seqs)
    [N, L] = size(seqs);
    order = matDim(E);
    indices = getIndeices1D(seqs, order);
    indices = reshape(indices, [L - order + 1, N]);
    logLikes = sum(log(E(indices)), 1);
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

