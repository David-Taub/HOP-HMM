function [accuricy, amounts] = learn(posSeqs, negSeqs, overlaps)
    % get the n'st most frequent overlap
    % load('/cs/stud/boogalla/  projects/CompGenetics/BaumWelch/peaks.mat');
    % posSeqs = seqs;
    % negSeqs = readSeq('NEnhancers.seq', L);
    % order = 1;
    for order = 2:5
        % order
        negCGTrain = true;
        negCGTest = false;
        M = size(overlaps, 2); %23
        accuricy = zeros(1, M + 1);
        amounts = zeros(1, M + 1);
        thresholds = zeros(1, M + 1);
        Es = zeros(2 * 4 ^ order, M + 1);
        freqDiffs = zeros(4 ^ order, M + 1);
        datasets = loadSeqs(posSeqs, negSeqs, overlaps, negCGTrain, negCGTest);
        % plotLettersFreq(datasets);
        freqRegression(datasets, order);
        % for overlapClass = 0:M
        %     j = overlapClass + 1;
        %     [MMResult, ~, freqDiffs(:, j), Es(:, j), thresholds(j), amounts(j)] = sampleAndLearnMulti(datasets{overlapClass + 1}, order);
        %     accuricy(j) = MMResult.ACC;
        %     MMResult
        % end

        % crossClassify(datasets, Es, thresholds, order);
        % diffHistPlot(freqDiffs, M + 1, shuffleNeg, overlaps);
        % crossLikelihood(datasets{1}.XTest(datasets{1}.YTest == 1, :), Es, datasets{1}.testOverlaps, order);
        % indicativeMotifsPlot(freqDiffs);
    end
end


function plotLettersFreq(datasets)
    figure
    p = [];
    M = length(datasets);
    for overlapClass = 1:M
        x = datasets{overlapClass}.XTrain(datasets{overlapClass}.YTrain == 1, :);
        p(1, overlapClass) = sum(x(:) == 1) / length(x(:));
        p(2, overlapClass) = sum(x(:) == 2) / length(x(:));
        p(3, overlapClass) = sum(x(:) == 3) / length(x(:));
        p(4, overlapClass) = sum(x(:) == 4) / length(x(:));
    end

    x = datasets{1}.XTrain(datasets{1}.YTrain == 2);
    p(1, M + 1) = sum(x(:) == 1) / length(x(:));
    p(2, M + 1) = sum(x(:) == 2) / length(x(:));
    p(3, M + 1) = sum(x(:) == 3) / length(x(:));
    p(4, M + 1) = sum(x(:) == 4) / length(x(:));
    [~, ord] = sort(var(p,[],1));
    
    hold on
    plot(p(1,ord))
    plot(p(2,ord))
    plot(p(3,ord))
    plot(p(4,ord))
    legend('A', 'C', 'G', 'T')
    addTissuesTicks(true, true, ord, false);

    
end 

function addTissuesTicks(withAll, withBackground, ord, withY)
    tissues = {'All', 'BAT', 'BMDM', 'BoneMarrow',...
           'CH12', 'Cerebellum', 'Cortex',...
           'E14', 'Heart-E14.5', 'Heart',...
           'Kidney', 'Limb-E14.5', 'Liver-E14.5',...
           'Liver', 'MEF', 'MEL', 'OlfactBulb',...
           'Placenta', 'SmIntestine', 'Spleen',...
           'Testis', 'Thymus', 'WholeBrain-E14.5', 'mESC', 'background'};
    if ~withAll
        tissues = tissues(2:end);
    end
    if ~withBackground
        tissues = tissues(1:end-1);
    end
    M = length(tissues);
    if ord == 0
        ord = 1:M;
    end
    ax = gca;
    ax.XTick = 1:M;
    ax.XTickLabel = tissues(ord);
    ax.XTickLabelRotation=45;
    if withY
        ax.YTick = 1:M;
        ax.YTickLabel = tissues(ord);
    end
end

function freqRegression(datasets, order)
    polynomialOrder = 1;
    M = size(datasets{1}.trainOverlaps, 2);
    Rsquare = zeros(1,M);
    % for tissue = 1 : M;
    for tissue = 1 : M;
        % N x 1
        trainOverlaps = log(1 + datasets{tissue + 1}.trainOverlaps(:, tissue));
        testOverlaps = log(1 + datasets{tissue + 1}.testOverlaps(:, tissue));
        trainOverlaps = [trainOverlaps; testOverlaps];
        % N x L
        trainPos = datasets{tissue + 1}.XTrain(datasets{tissue + 1}.YTrain == 1, :);
        % trainNeg = datasets{tissue + 1}.XTrain(datasets{tissue + 1}.YTrain == 2, :);
        testPos = datasets{tissue + 1}.XTrain(datasets{tissue + 1}.YTest == 1, :);
        trainPos = [trainPos; testPos];
        [NTrain, ~] = size(trainPos);
        [NTest, L] = size(testPos);


        % N x (L - order + 1)
        trainPosI = reshape(getIndeices1D(trainPos, order), [L - order + 1, NTrain]).';
        % trainNegI = reshape(getIndeices1D(trainNeg, order), [L - order + 1, NTrain]).';
        testPosI = reshape(getIndeices1D(testPos, order), [L - order + 1, NTest]).';

        % regressors
        k = 4 ^ order; %number of possible motifs with length 'order'
        % N x k
        trainPosHist = normr(histc(trainPosI, 1 : k, 2));
        % trainNegHist = normr(histc(trainNegI, 1 : k, 2));
        testPosHist = normr(histc(testPosI, 1 : k, 2));


        % % feature selection:
        % newK = 100 ;
        % newK = min(newK, k) ;
        % bestMotifs = selectBestMotifs(trainNegHist, trainPosHist, newK);
        % N x newK
        % trainPosHist = trainPosHist(:, bestMotifs);
        % testPosHist = testPosHist(:, bestMotifs);
        
        
        % % N x newK * polynomialOrder + 1
        % trainRegressor = zeros(NTrain, polynomialOrder * newK  + 1);
        % testRegressor = zeros(NTest, polynomialOrder * newK  + 1);

        % % make polynomial
        % for j = 1:polynomialOrder
        %     trainRegressor(:, newK  * (j-1) + 1 : newK  * j) = trainPosHist .^ j;
        %     testRegressor(:, newK  * (j-1) + 1 : newK  * j) = testPosHist .^ j;
        % end
        % % N x newK * polynomialOrder + 1
        % trainRegressor(:, end) = 1;
        % testRegressor(:, end) = 1;

        
        % linear regression
        % newK * polynomialOrder + 1 x 1

        model = 'linear';
        size(trainPosHist)
        if size(trainPosHist, 1) < size(trainPosHist, 2)
            continue;
        end
        stats = regstats(trainOverlaps, trainPosHist, model);

        Rsquare(tissue) = stats.rsquare
        % figure;
        % scatter(trainOverlaps, stats.yhat);

        matSize = [4 * ones(1, order), 1];
        N =4^order;

        vec = 1:N;

        a = zeros(N, length(matSize));
        for i = 1:length(matSize)
            a(:, i) = mod(vec, matSize(i));
            vec = ceil(vec / matSize(i));
        end

        % [a(:, 1), a(:, 2), a(:, 3), a(:, 4), a(:, 5)] = ind2sub(matSize, 1:N);
        for i = 1:N
            let{i} = int2nt(a(i, :));
        end
        % figure;
        % plot(stats.beta)
        % grid on;
        % ax = gca;
        % ax.XTick = 1:N;
        % ax.XTickLabel = let;
        % ax.XTickLabelRotation=45;
            
        % x = trainRegressor \ trainOverlaps;
        % % least square (not working)
        % x = (trainRegressor.' * trainRegressor) \ (trainRegressor.' * trainOverlaps);
        % alpha = 100;
        % x = (trainRegressor.' * trainRegressor + eye(size(trainRegressor, 2)) .* alpha) \ (trainRegressor.' * trainOverlaps);

        % [~, ord] = sortrows(testOverlaps);
        % testOverlaps = testOverlaps(ord, :);
        % estTestOverlaps = estTestOverlaps(ord, :);
        % subplot(1,2,1); imagesc(testOverlaps); colorbar;
        % subplot(1,2,2); imagesc(estTestOverlaps); colorbar;


    end
    hold on
    Rsquare(Rsquare<0) = 0
    plot(Rsquare)
    addTissuesTicks(false, false, 0, false)
end


% best motifs are the ones which are the most stably different from background
function bestMotifs = selectBestMotifs(trainNegHist, trainPosHist, N)
    diff = trainNegHist - trainPosHist;
    diffMean = abs(mean(diff, 1));
    diffVar = var(diff, [], 1);
    diffMeanSubMax = max(diffMean) - diffMean;
    [~,ord] = sort(diffMeanSubMax .^ 2 + diffVar .^ 2, 'ascend');
    bestMotifs = ord(1:N);
end
function crossLikelihood(posSeqs, Es, overlaps, order)
    [N, M] = size(overlaps);
    logLikes = zeros(N, M);
    for j = 1:M
        E = reshape(Es(:, j+1), [2, 4 .* ones(1, order)]);
        if order > 1
            posE = reshape(E(1, :), 4 .* ones(1, order));
        else
            posE = E(1, :);
        end
        logLikes(:, j) = getLogLikes(posE, posSeqs);
    end
    maxLogLikes = max(logLikes, [], 2);
    logLikes = bsxfun(@minus, logLikes, maxLogLikes);
    logLikes = exp(logLikes);

    % reorder
    [~, ord] = sortrows(overlaps);
    overlaps = overlaps(ord, :);
    logLikes = logLikes(ord, :);

    diffLogLikes = logLikes - (overlaps > 0);

    fig = figure();
    subplot(1,3,1); imagesc(diffLogLikes); colorbar;
    addTissuesTicks(false, false, 0, false);
    title('likelihood - overlaps');
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
    saveas(fig,sprintf('/cs/cbio/david/projects/CompGenetics/mm9Genome/graphs/likelihood/O%d_SH_SH.jpg', order));
end

% sample datasets, train and get test errors multiple times
% also get motifs frequency analyzed
function crossClassify(datasets, Es, thresholds, order)
    fprintf('order %d\n', order);
    M = length(datasets);
    crossClassMat = zeros(M, M);
    parfor j = 1:M
        E = reshape(Es(:, j), [2,ones(1,order) * 4]);
        threshold = thresholds(j);
        for i = 1:M
            [result, ~] = classify(E, datasets{i}.XTest, datasets{i}.YTest, threshold);
            crossClassMat(j, i) = result.ACC;
            fprintf('order %d %d %d %f \n', order, j, i, result.ACC);
            
            % plotHeatMap(crossClassMat(:,:), true, false, true, 'Error Map');
            % drawnow
        end
    end
    fig = plotHeatMap(crossClassMat, true, true, true, sprintf('Error Map: Models vs. Datasets (%d)', order));
    saveas(fig,sprintf('/cs/cbio/david/projects/CompGenetics/mm9Genome/graphs/errorMap/O%d_SH_SH.jpg', order));
    drawnow;
end

function f = plotHeatMap(tissueDat, reorder, dendro, isSquare, plotTitle)
    f = 0;
    M = size(tissueDat, 2);
    
    % reorder
    if reorder
        tree = linkage(tissueDat);
        leafOrd = optimalleaforder(tree, pdist(tissueDat));
        tissueDat = tissueDat(leafOrd, :);
        if isSquare
            tissueDat = tissueDat(:, leafOrd);
        end
        if dendro
            figure;
            dendrogram(tree,'Reorder',leafOrd)
            addTissuesTicks(true, false, leafOrd, false);
            title(plotTitle);
            f = figure();
        end
    end

    if isSquare
        imagesc(tissueDat, [0,1]);colorbar;
    else
        D = squareform(pdist(tissueDat));
        imagesc(D,[0,1]);colorbar;
    end
    title(plotTitle);
    addTissuesTicks(true, false, leafOrd, true);

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

function diffHist = freqFinder(posSeqs, negSeqs, order)
    indicesP = getIndeices1D(posSeqs, order);
    indicesN = getIndeices1D(negSeqs, order);
    posHist = histc(indicesP, 1 : 4 ^ order);
    negHist = histc(indicesN, 1 : 4 ^ order);
    diffHist = log(posHist/negHist);
    diffHist = diffHist / size(posSeqs, 1);
    % plot(sort(diffHist))
    % hold on
end
function indicativeMotifsPlot(diffHist)
    [~, sortedDiff] = sort(diffHist, 1);
    Ps = 1:floor(length(diffHist) / 2);
    vals = zeros(1, length(Ps));
    for i = 1:length(Ps)
        curS = sortedDiff;
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
        
        saveas(fig,filepath);
    end
end

function out = getReverseComplement(seqs)
    out = fliplr(5 - seqs);
end

% negCGTest may be true only when negCGTrain is true
function datasets = loadSeqs(posSeqs, negSeqs, overlaps, negCGTrain, negCGTest)
    M = size(overlaps, 2);
    N = min(size(negSeqs, 1), size(posSeqs, 1));
    datasets = cell(M+1, 1);
    trainTestRate = 0.8;
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
    datasets{1}.trainOverlaps = overlaps(posOrder(1:trainLabLength), :);
    datasets{1}.testOverlaps = overlaps(posOrder(trainLabLength+1:N), :);
    
    for i = 1 : M
        trainPos = datasets{1}.trainOverlaps(:, i) > 0;
        testPos = datasets{1}.testOverlaps(:, i) > 0;
        datasets{i + 1}.XTrain = datasets{1}.XTrain([trainPos; trainPos], :);
        datasets{i + 1}.XTest = datasets{1}.XTest([testPos; testPos], :);
        datasets{i + 1}.YTrain = datasets{1}.YTrain([trainPos; trainPos]);
        datasets{i + 1}.YTest = datasets{1}.YTest([testPos; testPos]);
        datasets{i + 1}.trainOverlaps = datasets{1}.trainOverlaps(trainPos, :);
        datasets{i + 1}.testOverlaps = datasets{1}.testOverlaps(testPos, :);
    end

end

function [testErrMM, testErrHMM, freqDiff, E] = learnData(XTrain, YTrain, XTest, YTest, order)
    E = trainMarkov(XTrain, YTrain, order);
    thresholds = 0.8 : 0.005 : 1.2;
    [~, threshold] = classify(E, XTrain, YTrain, thresholds);
    [testErrMM, ~] = classify(E, XTest, YTest, threshold);
    freqDiff = freqFinder(XTrain(YTrain == 1, :), XTrain(YTrain == 2, :), order);
    testErrHMM = 0;
    % testErrMM = 0;

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
% N - number of sequences to calculate the posterior with
% out - N x L
function out = getPosterior(seqs, startT, T, E)
    m = 2;
    [N, L] = size(seqs);
    posterior = zeros(m, L, N);
    [alpha, scale] = forwardAlg(seqs, startT, T, E);
    beta = backwardAlg(seqs, startT, T, E, scale);
    % S x m x L
    posterior = alpha .* beta;
    posterior = bsxfun(@times, posterior, 1 ./ sum(posterior, 2));
    % return posterior of the positive state, 
    out(:,:) = posterior(:, 1, :);
end

function [startT, T] = createHmmParams(neg2pos, pos2neg)
    T = [1 - pos2neg, pos2neg; neg2pos, 1 - neg2pos];
    startT = [0.5; 0.5];
end

function E = trainMarkov(X, Y, order)
    E = [];
    % for i = [unique(Y)]
    for i = 1:max(Y, [], 1)
        seqs = [X(Y == i, :); getReverseComplement(X(Y == i, :))];
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
    else
        threshold = thresholds;
    end
    err = getLose(ratioPos.', ratioNeg.', threshold);
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


function anaFreq(posSeqs, negSeqs, order)
    matSize = [4 * ones(1, order), 1];
    indicesP = getIndeices1D(posSeqs, order);
    indicesN = getIndeices1D(negSeqs, order);
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

