
function [accuricy, amounts] = learn(posSeqs, negSeqs, overlaps)
    % get the n'st most frequent overlap
    % load('/cs/stud/boogalla/  projects/CompGenetics/BaumWelch/peaks.mat');
    % posSeqs = seqs;
    % negSeqs = readSeq('NEnhancers.seq', L);
    tissues = {'all', 'BAT', 'BMDM', 'BoneMarrow',...
               'CH12', 'Cerebellum', 'Cortex',...
               'E14', 'Heart-E14.5', 'Heart',...
               'Kidney', 'Limb-E14.5', 'Liver-E14.5',...
               'Liver', 'MEF', 'MEL', 'OlfactBulb',...
               'Placenta', 'SmIntestine', 'Spleen',...
               'Testis', 'Thymus', 'WholeBrain-E14.5', 'mESC'};
    order = 6;
    negCGTrain = false;
    negCGTest = false;
    reorderHeatMap = true;
    M = size(overlaps, 2); %23
    accuricy = zeros(1, M+1);
    amounts = zeros(1, M+1);
    thresholds = zeros(1, M+1);
    Es = zeros(2 * 4 ^ order, M + 1);
    freqDiffs = zeros(4 ^ order, M + 1);
    datasets = loadSeqs(posSeqs, negSeqs, overlaps, negCGTrain, negCGTest);
    for overlapClass = 0:0%M
        j = overlapClass + 1;
        [MMResult, ~, freqDiffs(:, j), Es(:, j), thresholds(j), amounts(j)] = sampleAndLearnMulti(datasets{overlapClass + 1}, order);
        accuricy(j) = MMResult.ACC;
    end

    % crossClassify(datasets, Es, thresholds, order);
    % diffHistPlot(freqDiffs, M + 1, shuffleNeg, reorderHeatMap);
    % crossLikelihood(datasets, Es, overlaps, order);
    % indicativeMotifsPlot(freqDiffs);
end



function crossLikelihood(posSeqs, Es, overlaps, order)
    [N, M] = size(overlaps);
    logLikes = zeros(N, M);
    for j = 1:M
        E = reshape(Es(:, j+1), [2, 4 .* ones(1, order)]);
        posE = reshape(E(1, :), 4 .* ones(1, order));
        logLikes(:, j) = getLogLikes(posE, posSeqs);
    end
    maxLogLikes = max(logLikes, [], 2);
    logLikes = bsxfun(@minus, logLikes, maxLogLikes);
    logLikes = exp(logLikes);

    % reorder
    [~, ord] = sortrows(overlaps);
    overlaps = overlaps(ord, :);
    logLikes = logLikes(ord, :);

    diffLogLikes = logLikes - overlaps;

    figure
    subplot(1,3,1); imagesc(diffLogLikes); colorbar;
    tissues = {'BAT', 'BMDM', 'BoneMarrow',...
               'CH12', 'Cerebellum', 'Cortex',...
               'E14', 'Heart-E14.5', 'Heart',...
               'Kidney', 'Limb-E14.5', 'Liver-E14.5',...
               'Liver', 'MEF', 'MEL', 'OlfactBulb',...
               'Placenta', 'SmIntestine', 'Spleen',...
               'Testis', 'Thymus', 'WholeBrain-E14.5', 'mESC'};
    ax = gca;
    ax.XTick = 1:M;
    ax.XTickLabel = tissues;
    ax.XTickLabelRotation=45;
    title('diff');
    subplot(1,3,2); imagesc(overlaps); colorbar;
    ax = gca;
    ax.XTick = 1:M;
    ax.XTickLabel = tissues;
    ax.XTickLabelRotation=45;
    title('overlaps');
    subplot(1,3,3); imagesc(logLikes); colorbar;
    ax = gca;
    ax.XTick = 1:M;
    ax.XTickLabel = tissues;
    ax.XTickLabelRotation=45;
    title('likelihood');
    
end

% sample datasets, train and get test errors multiple times
% also get motifs frequency analyzed
function crossClassify(datasets, Es, thresholds, order)
    M = length(datasets);
    crossClassMat = zeros(M, M, repeats);
    for j = 1:M
        E = reshape(Es(:, j + 1), [2,ones(1,order) * 4]);
        threshold = thresholds(j);
        for i = 1:M
            [result, ~] = classify(E, datasets{j}.XTest, datasets{j}.YTest, threshold);
            crossClassMat(j, i) = result.ACC;
            
            plotHeatMap(crossClassMat(:,:), true, false, true, 'Error Map');
            drawnow
        end
    end
    plotHeatMap(crossClassMat, true, true, true, 'Error Map: All Tissue Models vs. Tissue Datasets');
end

function f = plotHeatMap(tissueDat, reorder, dendro, isSquare, plotTitle)
    f = 0;
    tissues = {'all', 'BAT', 'BMDM', 'BoneMarrow',...
               'CH12', 'Cerebellum', 'Cortex',...
               'E14', 'Heart-E14.5', 'Heart',...
               'Kidney', 'Limb-E14.5', 'Liver-E14.5',...
               'Liver', 'MEF', 'MEL', 'OlfactBulb',...
               'Placenta', 'SmIntestine', 'Spleen',...
               'Testis', 'Thymus', 'WholeBrain-E14.5', 'mESC'};
    M = length(tissues);
    
    % reorder
    if reorder
        tree = linkage(tissueDat);
        leafOrd = optimalleaforder(tree, pdist(tissueDat));
        tissues = tissues(leafOrd);
        tissueDat = tissueDat(leafOrd, :);
        if isSquare
            tissueDat = tissueDat(:, leafOrd);
        end
        if dendro
            figure;
            dendrogram(tree,'Reorder',leafOrd)
            ax = gca;
            ax.XTick = 1:M;
            ax.XTickLabel = tissues;
            ax.XTickLabelRotation=45;
            title(plotTitle);
            f = figure();
        end
    end

    if isSquare
        imagesc(tissueDat);colorbar;
    else
        D = squareform(pdist(tissueDat));
        imagesc(D);colorbar;
    end
    title(plotTitle);
    ax = gca;
    ax.XLim = [0.5 M+0.5];
    ax.YLim = [0.5 M+0.5];
    ax.XTick = 1:M;
    ax.YTick = 1:M;
    ax.YTickLabel = tissues;
    ax.XTickLabel = tissues;
    ax.XTickLabelRotation=45;
    title(plotTitle);

end 
% sample datasets, train and get test errors multiple times
% also get motifs frequency analyzed
function [MMErr, HMMErr, freqDiff, E, threshold, amount] = sampleAndLearnMulti(dataset, order)


    thresholdRange = 0.8 : 0.005 : 1.2;
    amount = sum([dataset.YTrain == 1;dataset.YTest == 1], 1);
    % train
    E = trainMarkov(dataset.XTrain, dataset.YTrain, order);
    
    % classify train and test datasets
    [~, threshold] = classify(E, dataset.XTrain, dataset.YTrain, thresholdRange);
    [MMErr, ~] = classify(E, dataset.XTest, dataset.YTest, threshold);

    % post analyze
    freqDiff = freqFinder(dataset.XTrain(dataset.YTrain == 1, :), dataset.XTrain(dataset.YTrain == 2, :), order);
    % [testErrMMs(i), testErrHMMs(i), freqDiffs(:, i), E] = learnData(dataset.XTrain, YTrain, dataset.XTest, dataset.YTest, order);
    % freqDiffs(:, :, i) = spread(reshape(E(1, :), ones(1, order) * 4));
    HMMErr = 0;
    E = E(:);
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
    % diffHist = posHist - negHist;
    % diffHist(diffHist < 0) = 0;
    diffHist = diffHist / size(seqsPos, 1);
    % plot(sort(diffHist))
    % hold on
end
function indicativeMotifsPlot(diffHist)
    [~, s] = sort(diffHist, 1);
    size(s)
    N = length(diffHist);
    Ps = 1:floor(N / 2);
    vals = zeros(1, length(Ps));
    for i = 1:length(Ps)
        curS = s;
        curS(Ps(i):end - Ps(i), :) = [];
        vals(i) = length(unique(curS(:)));
    end
    figure 
    plot(Ps, vals);
end
% Es = 4 ^ order x M
function diffHistPlot(diffHist, M, shuffleNeg, reorderHeatMap)
    for p = [1,2,3,5,10,20,200]
        [~, s] = sort(diffHist, 1);
        s(p:end - p, :) = [];
        diffHistU = diffHist(unique(s(:)), :);

        fig = plotHeatMap(diffHistU.', reorderHeatMap, true, false, sprintf('Emission Diff Between Tissues %d %d', p, size(diffHistU,1)));

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
        
        saveas(fig,filepath);%close all;
    end
end

% negCGTest may be true only when negCGTrain is true
function datasets = loadSeqs(posSeqs, negSeqs, overlaps, negCGTrain, negCGTest)
    M = size(overlaps, 2);
    N = min(size(negSeqs, 1), size(posSeqs, 1));
    datasets = cell(M+1, 1);
    trainTestRate = 0.9;
    trainLabLength = ceil(N * trainTestRate);

    % shuffle pos
    posOrder = randperm(size(posSeqs, 1), N);
    posSeqs = posSeqs(posOrder, :);
    

    % shuffle neg
    % Note: the negative are given sorted by their ACGT content, so taking only the top means 
    % the negative seqs are most alike the ACGT content of enhancers
    if negCGTrain
        if negCGTest
            negSeqs = negSeqs(randperm(N, N), :);
        else
            mostCG = negSeqs(1:trainLabLength, :);
            leastCG = negSeqs(trainLabLength + 1:end, :);
            negSeqs = [mostCG(randperm(trainLabLength, trainLabLength), :);...
                       leastCG(randperm(size(leastCG, 1), N - trainLabLength), :) ];
        end
    else
        negSeqs = negSeqs(randperm(size(negSeqs, 1), N), :);
    end
    
    % build 'all' (1) dataset
    datasets{1}.XTrain = [posSeqs(1:trainLabLength, :); negSeqs(1:trainLabLength, :)];
    datasets{1}.XTest  = [posSeqs(trainLabLength + 1: N, :); negSeqs(trainLabLength + 1:N, :)];
    datasets{1}.YTrain = [ones(trainLabLength, 1); ones(trainLabLength, 1) .* 2];
    datasets{1}.YTest  = [ones(N - trainLabLength,1); ones(N - trainLabLength,1) .* 2];
    datasets{1}.overlapsTrain = overlaps(posOrder(1:trainLabLength), :);
    datasets{1}.overlapsTest = overlaps(posOrder(trainLabLength+1:N), :);
    
    for i = 1 : M
        trainPos = datasets{1}.overlapsTrain(:, i) == 1;
        testPos = datasets{1}.overlapsTest(:, i) == 1;
        datasets{i + 1}.XTrain = datasets{1}.XTrain([trainPos; trainPos], :);
        datasets{i + 1}.XTest = datasets{1}.XTest([testPos; testPos], :);
        datasets{i + 1}.YTrain = datasets{1}.YTrain([trainPos; trainPos]);
        datasets{i + 1}.YTest = datasets{1}.YTest([testPos; testPos]);
        datasets{i + 1}.overlapsTrain = datasets{1}.overlapsTrain(trainPos, :);
        datasets{i + 1}.overlapsTest = datasets{1}.overlapsTest(testPos, :);
    end

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
function results = getLose(high, low, threshold)
    
    TP = sum(high > threshold, 1);
    FP = sum(low  > threshold, 1);
    FN = sum(high < threshold, 1);
    TN = sum(low  < threshold, 1);
    
    results.MEA = length(high);
    results.ACC = (TN + TP) / (FN + TN + FP + TN);
    results.TPR = (TP) / (TP + FN); % sensitivity \ recall
    results.FNR = (FN) / (TP + FN); % miss rate
    results.TNR = (TN) / (TN + FP); % specificity
    results.FPR = (FP) / (TN + FP); % fall out
    results.MCC = (TP * TN - FP * FN) / sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN));
    
    % AUC ROC
    thresholdRange = min([high;low]) : 0.001 : max([high;low]);
    L = length(thresholdRange);
    % 1 x L
    TPs = sum(bsxfun(@gt, repmat(high, [1, L]), thresholdRange), 1);
    FNs = sum(bsxfun(@lt, repmat(high, [1, L]), thresholdRange), 1);
    FPs = sum(bsxfun(@gt, repmat(low, [1, L]), thresholdRange), 1);
    TNs = sum(bsxfun(@lt, repmat(low, [1, L]), thresholdRange), 1);
    TPRs = (TPs) ./ (TPs + FNs);
    FPRs = (FPs) ./ (TNs + FPs);
    results.AUC = trapz(FPRs(end:-1:1), TPRs(end:-1:1));
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
    ratioPos = getLikeRatio(E, X(Y == 1, :)); %high
    ratioNeg = getLikeRatio(E, X(Y == 2, :)); %low
    if length(thresholds) > 1
        [~, threshold] = findThreshold(ratioNeg.', ratioPos.', thresholds);
        err = getLose(ratioPos.', ratioNeg.', threshold);
    else
        threshold = thresholds;
        err = getLose(ratioPos.', ratioNeg.', threshold);
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

