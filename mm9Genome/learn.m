function [MMmean, amounts] = learn(posSeqs, negSeqs, overlaps)
    % get the n'st most frequent overlap
    % load('/cs/stud/boogalla/  projects/CompGenetics/BaumWelch/peaks.mat');
    % posSeqs = seqs;
    % negSeqs = readSeq('NEnhancers.seq', L);
    order = 6;
    shuffleNeg = true;
    reorderHeatMap = true;
    M = size(overlaps, 2); %23
    MMmean = zeros(1, M+1);
    amounts = zeros(1, M+1);
    Es = zeros(2 * 4 ^ order, M + 1);
    freqDiffs = zeros(4 ^ order, M + 1);
    for overlapClass = 0:M
        j = overlapClass + 1;
        % [MMmean(overlapClass+1), ~, Es(:, :, overlapClass + 1), amounts(overlapClass+1)] = sampleAndLearnMulti(posSeqs, negSeqs, overlaps, overlapClass, order);
        [MMmean(j), ~, freqDiffs(:, j), Es(:, j), amounts(j)] = sampleAndLearnMulti(posSeqs, negSeqs, overlaps, overlapClass, order, shuffleNeg);
    end
    % diffHistPlot(freqDiffs, M + 1, shuffleNeg, reorderHeatMap);
end


% sample datasets, train and get test errors multiple times
% also get motifs frequency analyzed
function [MMmean, HMMmean, freqDiff, E, amount] = sampleAndLearnMulti(posSeqs, negSeqs, overlaps, overlapClass, order, shuffleNeg)

    repeats = 2;

    testErrMMs  = zeros(1, repeats);
    % testErrHMMs = zeros(1, repeats);
    freqDiffs = zeros(4 ^ order, repeats);
    Es = zeros(2 * 4 ^ order, repeats);
    for i = 1:repeats
        [XTrain, XTest, YTrain, YTest] = loadSeqs(posSeqs, negSeqs, overlaps, overlapClass, shuffleNeg);
        amount = sum([YTrain == 1;YTest == 1], 1);
        % train
        E = trainMarkov(XTrain, YTrain, order);
        thresholds = 0.8 : 0.005 : 1.2;
        
        % classify train and test datasets
        [~, threshold] = classify(E, XTrain, YTrain, thresholds);
        [testErrMMs(i), ~] = classify(E, XTest, YTest, threshold);

        % post analyze
        freqDiffs(:, i) = freqFinder(XTrain(YTrain == 1, :), XTrain(YTrain == 2, :), order);
        % [testErrMMs(i), testErrHMMs(i), freqDiffs(:, i), E] = learnData(XTrain, YTrain, XTest, YTest, order);
        Es(:, i) = E(:);
        % freqDiffs(:, :, i) = spread(reshape(E(1, :), ones(1, order) * 4));
    end
    MMmean = mean(testErrMMs, 2);
    [order, MMmean, amount]
    % HMMmean = mean(testErrHMMs, 2);
    HMMmean = 0;
    freqDiff = mean(freqDiffs, 2);
    E = mean(Es, 2);
end

% E - 4 x 4 x ... x 4 ('order' times)
% outE - 4 x n  where n = 4 ^ (order -1)
% function outE = spread(E)
%     order = matDim(E);
%     outE = shiftdim(E, order-1);
%     outE = outE(:,:);
% end

% function plotEs(Es, M)
%     tissues = {'BAT', 'BMDM', 'BoneMarrow',...
%                'CH12', 'Cerebellum', 'Cortex',...
%                'E14', 'Heart-E14.5', 'Heart',...
%                'Kidney', 'Limb-E14.5', 'Liver-E14.5',...
%                'Liver', 'MEF', 'MEL', 'OlfactBulb',...
%                'Placenta', 'SmIntestine', 'Spleen',...
%                'Testis', 'Thymus', 'WholeBrain-E14.5', 'mESC'};

%     Es = mean(Es, 2);
%     EsRep = repmat(Es, [1, M, 1]);
%     %todo: remove cells without any instances get top X% of different cells
%     D = median((EsRep - permute(EsRep, [2,1,3])) .^ 2, 3);
%     figure;
%     imagesc(D);colorbar;
%     title('Emission Diff Between Tissues');
%     % set(gca,'YLim',[0 M],'YTick',1:12,'YTickLabel',months)
%     ax = gca;
%     ax.XLim = [0 M];
%     ax.YLim = [0 M];
%     ax.XTick = 1:M;
%     ax.YTick = 1:M;
%     ax.YTickLabel = tissues;
%     ax.XTickLabel = tissues;
%     ax.XTickLabelRotation=45;
%     % set(gca,'YLim',[0 M],'YTick',1:M, 'XLim',[0 M],'XTick',1:M,...
%     %         'YTickLabel', tissues, 'XTickLabel', tissues);
% end

function diffHist = freqFinder(seqsPos, seqsNeg, order)
    indicesP = getIndeices1D(seqsPos, order);
    indicesN = getIndeices1D(seqsNeg, order);
    posHist = histc(indicesP, 1 : 4 ^ order);
    negHist = histc(indicesN, 1 : 4 ^ order);
    diffHist = log(posHist/negHist);
    diffHist = posHist - negHist;
    % diffHist(diffHist < 0) = 0;
    diffHist = diffHist / size(seqsPos, 1);
    plot(sort(diffHist))
    hold on
end

% Es = 4 ^ order x M
function diffHistPlot(diffHist, M, shuffleNeg, reorderHeatMap)
    for p = [1,2,3,5,10,20,200]
        tissues = {'all', 'BAT', 'BMDM', 'BoneMarrow',...
                   'CH12', 'Cerebellum', 'Cortex',...
                   'E14', 'Heart-E14.5', 'Heart',...
                   'Kidney', 'Limb-E14.5', 'Liver-E14.5',...
                   'Liver', 'MEF', 'MEL', 'OlfactBulb',...
                   'Placenta', 'SmIntestine', 'Spleen',...
                   'Testis', 'Thymus', 'WholeBrain-E14.5', 'mESC'};
        [~, s] = sort(diffHist, 1);
        s(p:end - p, :) = [];
        diffHistU = diffHist(unique(s(:)), :);

        tree = linkage(diffHistU.');
        
        % reorder
        if reorderHeatMap
            leafOrd = optimalleaforder(tree, pdist(diffHistU.'));
            tissues = tissues(leafOrd);
            diffHistU = diffHistU(:, leafOrd);
        end
        dendrogram(tree);
        distMat = squareform(pdist(diffHistU.') .^ 2);

        f = figure;
        imagesc(distMat);colorbar;
        title(sprintf('Emission Diff Between Tissues %d %d', p, size(diffHistU,1)));
        if reorderHeatMap
            if shuffleNeg
                filepath = sprintf('/a/store-05/z/cbio/david/projects/CompGenetics/mm9Genome/graphs/motifs/reordered/Motifs_%d_S.jpg', p);
            else
                filepath = sprintf('/a/store-05/z/cbio/david/projects/CompGenetics/mm9Genome/graphs/motifs/reordered/Motifs_%d.jpg', p);
            end
        else
            if shuffleNeg
                filepath = sprintf('/a/store-05/z/cbio/david/projects/CompGenetics/mm9Genome/graphs/motifs/Motifs_%d_S.jpg', p);
            else
                filepath = sprintf('/a/store-05/z/cbio/david/projects/CompGenetics/mm9Genome/graphs/motifs/Motifs_%d.jpg', p);
            end
        end
        % set(gca,'YLim',[0 M],'YTick',1:12,'YTickLabel',months)
        ax = gca;
        ax.XLim = [0 M];
        ax.YLim = [0 M];
        ax.XTick = 1:M;
        ax.YTick = 1:M;
        ax.YTickLabel = tissues;
        ax.XTickLabel = tissues;
        ax.XTickLabelRotation=45;
        saveas(f,filepath);%close all;
    end
end

% L - sequence lengths
% n - number of classes to get for the positive sequences, where class
% means unique overlap between tissues, and if n is 1:3 then we take 
% the sequences of the three most frequent class
function [XTrain, XTest, YTrain, YTest] = loadSeqs(posSeqs, negSeqs, overlaps, overlapClass, shuffleNeg)
    if overlapClass > 0
        posSeqs = posSeqs(overlaps(:, overlapClass) == 1, :);
    end

    % assumes more negative than positives
    N = min(size(negSeqs, 1), size(posSeqs, 1));
    trainTestRate = 0.9;
    trainLabLength = ceil(N * trainTestRate);

    % shuffle
    posSeqs = posSeqs(randperm(size(posSeqs, 1), N), :);
    if shuffleNeg
        negSeqs = negSeqs(randperm(size(negSeqs, 1), N), :);
    end

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

function [testErrMM, testErrHMM, freqDiff, E] = learnData(XTrain, YTrain, XTest, YTest, order)
    E = trainMarkov(XTrain, YTrain, order);
    thresholds = 0.8 : 0.005 : 1.2;
    [~, threshold] = classify(E, XTrain, YTrain, thresholds);
    [testErrMM, ~] = classify(E, XTest, YTest, threshold);
    [order, threshold, testErrMM]
    size(XTrain, 1) / 2
    freqDiff = freqFinder(XTrain(YTrain == 1, :), XTrain(YTrain == 2, :), order);
    testErrHMM = 0;
    % testErrMM = 0;

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
    ratioPos = getLikeRatio(E, X(Y == 1, :));
    ratioNeg = getLikeRatio(E, X(Y == 2, :));
        if length(thresholds) > 1
        [err, threshold] = findThreshold(ratioNeg.', ratioPos.', thresholds);
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
    posHist = histc(indicesP, 1 : 4 ^ order);
    negHist = histc(indicesN, 1 : 4 ^ order);
    diffHist = posHist - negHist;
    N = 40;
    [sortedX,sortingIndices] = sort(diffHist,'descend');
    % maxValues = sortedX(1:N);
    maxValueIndices = sortingIndices(1:N);
    maxValues = diffHist(maxValueIndices);
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
    elseif order == 3
        [a(:, 1), a(:, 2), a(:, 3)] = ...
        ind2sub(matSize, maxValueIndices);
    elseif order == 2
        [a(:, 1), a(:, 2)] = ...
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

