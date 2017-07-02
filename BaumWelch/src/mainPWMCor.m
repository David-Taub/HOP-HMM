    % checks correlation between the existence of p300 peak and its TF binding sites
% presence. the binding sights presence are evaluated with PWMs.

% 10  TxEnh5 Transcribed 5' preferential and Enh
% 11  TxEnh3 Transcribed 3' preferential and Enh
% 12  TxEnhW  Transcribed and Weak Enhancer
% 13  EnhA1   Active Enhancer 1
% 14  EnhA2   Active Enhancer 2
% 15  EnhAF   Active Enhancer Flank
% 16  EnhW1   Weak Enhancer 1
% 17  EnhW2   Weak Enhancer 2
% 18  EnhAc   Primary H3K27ac possible Enhancer

% usage:
% cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src
% mergedPeaksMin = load(fullfile('data', 'peaks', 'mergedPeaksMinimized.mat'));
% delete(fullfile('data', 'precomputation', 'pcPWMp.mat'));
% mergedPeaksMin = load('data/RoadmapEnhancers.mat');


% mergedPeaksMin = load('data/peaks/roadmap/mergedPeaksMinimized.mat');
% mainPWMCor(mergedPeaksMin);
function [maxPeaks, overlaps] = mainPWMCor(mergedPeaksMin)
    close all;
    tic
    % N = size(mergedPeaksMin.seqs, 1);
    fprintf('rearranging data\n')
    % [overlaps, Xs] = genData(mergedPeaksMin);
    [overlaps, Xs, maxPeaks] = genData(mergedPeaksMin);
    % [N, L] = size(Xs);
    % N x k x L
    size(maxPeaks)
    k = size(maxPeaks, 2); %205
    r = size(overlaps, 2); %19
    % N x k
    % maxPeaks = max(pcPWMp, [], 3);
    % maxPeaks = pcPWMp;
    % maxPeaks = maxPeaks + nMax(pcPWMp, 2);
    % assert(not(any(isnan(maxPeaks(:)))))
    % maxPeaks = maxPeaks + mean(pcPWMp, 3);

    % maxPeaks = maxPeaks + nMax(pcPWMp, 3);

    % [meansDiff, ~] = getMeanDiff(maxPeaks, overlaps);
    fprintf('separation tests\n')

    % load(fullfile('data', 'outMainPWMCor.mat'));
    ranks = getRanks(maxPeaks, overlaps, r, k);
    outputPath = fullfile('data', 'outMainPWMCor.mat');
    fprintf('Saving processed data to file %s\n', outputPath)
    save(outputPath, 'ranks', 'maxPeaks', 'overlaps');
    showAllSep(ranks);
    fprintf('showing results 1\n')
    showData(maxPeaks, ranks, overlaps);
    showBestSep2(ranks, maxPeaks, overlaps, r);
    fprintf('showing results 2\n')
    showBestSep(ranks, maxPeaks, overlaps, 1);
    showBestSep(ranks, maxPeaks, overlaps, 2);
    showBestSep(ranks, maxPeaks, overlaps, 3);
    toc
end

% maxPeaks - N x k x M
% overlaps - N x r
% ranks - r x k x M
function [ranks] = getRanks(maxPeaks, overlaps, r, k)
    ranks = zeros(r, k, 2);
    i = 1
    % 1 x k
    for j = 1:k
        peaksIndicatorTFPos = maxPeaks(overlaps(:, i) > 0, j);
        peaksIndicatorTFNeg = maxPeaks(overlaps(:, i) == 0, j);

        if mean(peaksIndicatorTFPos) > mean(peaksIndicatorTFNeg);
            ranks(i, j ,2) = getAucRoc(peaksIndicatorTFPos, peaksIndicatorTFNeg, false);
        else
            ranks(i, j, 2) = getAucRoc(peaksIndicatorTFNeg, peaksIndicatorTFPos, false);
        end

        [~, ranks(i, j, 1)] = kstest2(peaksIndicatorTFPos, peaksIndicatorTFNeg);
        % ranks(i, j, 1) = 1.0 - ranks(i, j, 1);

        % [ranks(i, j), ~] = ranksum(peaksIndicatorTFPos, peaksIndicatorTFNeg);
        % ranks(i, j) = 1 - ranks(i, j);
        % ranks(i, j) = matUtils.emdTest(peaksIndicatorTFPos, peaksIndicatorTFNeg);
        fprintf('%d / %d. %d / %d. %.2f\n', i, r, j, k, ranks(i, j));
    end
    fprintf('\n')
end

function auc = getAucRoc(pos, neg, shouldPlot)
    scores = [pos;neg];
    labels = [ones(length(pos), 1); zeros(length(neg), 1)];
    [X, Y, ~, auc] = perfcurve(labels, scores, 1);
    if shouldPlot
        plot(X,Y)
        xlabel('False positive rate')
        ylabel('True positive rate')
    end
end

% maxPeaks - N x k
% ranks - r x k
function showBestSep(ranks, maxPeaks, overlaps, r)
    [~,~,names] = BaumWelchPWM.PWMs();
    i = 1;
    % for i =1:r
        % 1 x k
    for j = 1:r
        [bestRank, tissueIndicator] = max(ranks(1, :, 2), [], 2);
        ranks(i, tissueIndicator, 2) = -inf;
    end
    peaksIndicatorTFPos = maxPeaks(overlaps(:, i) > 0, tissueIndicator);
    peaksIndicatorTFNeg = maxPeaks(overlaps(:, i) == 0, tissueIndicator);
    % subplot(4,5,i);
    figure;
    subplot(1,2,1);
    h = histogram(peaksIndicatorTFPos, 50, 'Normalization', 'probability');
    hold on;
    histogram(peaksIndicatorTFNeg, h.BinEdges, 'Normalization', 'probability');
    title(['PWM ', names{tissueIndicator}, ' LogLikes. Rate: ', num2str(bestRank)]);
    legend('Enhancers from tissue', 'Enhancers from all other tissues')
    subplot(1,2,2);
    if mean(peaksIndicatorTFPos) > mean(peaksIndicatorTFNeg);
        getAucRoc(peaksIndicatorTFPos, peaksIndicatorTFNeg, true);
    else
        getAucRoc(peaksIndicatorTFNeg, peaksIndicatorTFPos, true);
    end
    % end
end

% maxPeaks - N x k
% ranks - r x k
function showBestSep2(ranks, maxPeaks, overlaps, r)
    mask = ranks(1, :, 2) > 0.58;
    maxPeaks2 = max(maxPeaks(:, mask), [], 2);
    peaksIndicatorTFPos = maxPeaks2(overlaps(:, 1) > 0);
    peaksIndicatorTFNeg = maxPeaks2(overlaps(:, 1) == 0);
    figure;
    subplot(1,2,2);
    auc = getAucRoc(peaksIndicatorTFPos, peaksIndicatorTFNeg, true);
    subplot(1,2,1);
    h = histogram(peaksIndicatorTFPos, 50, 'Normalization', 'probability');
    hold on;
    histogram(peaksIndicatorTFNeg, h.BinEdges, 'Normalization', 'probability');
    title(['PWM (max of best) LogLikes. Rate: ', num2str(auc)]);
    legend('Enhancers from tissue', 'Enhancers from all other tissues')


    % end
end

% ranks - r x k
function showAllSep(ranks)
    figure;
    scatter(-log(ranks(1, :, 1)), (ranks(1, :, 2)));
    xlabel('-log of KS2 p-value (higher is better)');
    ylabel('AUC ROC (higher is better)');
    title(['PWM separation rank']);
    % end
end

% maxPeaks - N x k
function [meansDiff, tissueIndicator] = getMeanDiff(maxPeaks, overlaps, r, k)
    meansDiff = zeros(r, k);
    tissueIndicator = zeros(r, 1);
    for i =1:r
        % 1 x k
        meansDiff(i, :) = getAbsDiffMean(maxPeaks, overlaps, i);

        [~, tissueIndicator(i)] = max(meansDiff(i, :), [], 2);
        % peaksIndicatorTFPos = maxPeaks(overlaps(:, i) == 1, tissueIndicator(i));
        % peaksIndicatorTFNeg = maxPeaks(overlaps(:, i) == 0, tissueIndicator(i));

        % figure
        % histogram(peaksIndicatorTFPos, N / 10, 'Normalization', 'probability');
        % hold on;
        % histogram(peaksIndicatorTFNeg, N / 10, 'Normalization', 'probability');
        % title('peaks height of best motif in all enhancers');
        % drawnow;
    end
end

function showData(maxPeaks, ranks, overlaps)
    figure;
    [~, motifInd] = sort(ranks(1, :, 1), 2);
    % [maxPeaks, motifInd] = matUtils.clustMatRows(maxPeaks');
    maxPeaks = maxPeaks(:, motifInd);
    % [maxPeaks, seqInd] = matUtils.clustMatRows(maxPeaks');
    subplot(1,3,1);imagesc(maxPeaks); colorbar;
    xlabel('TF'); ylabel('Sequences');
    title('Max PSSM id difference sequence');


    % [ranks, tissueInd]  = matUtils.clustMatRows(ranks);
    ranks = shiftdim(ranks(1, motifInd, :), 1);

    subplot(1,3,2);imagesc(ranks(:,:,1)); colorbar;
    xlabel('TF'); ylabel('Tissue');
    title('AUC ROC of max PSSM of TF PWM (higher classifies better)')


    % overlaps = overlaps(seqInd, :);
    % overlaps = overlaps(:, tissueInd);
    subplot(1,3,3);imagesc(overlaps); colorbar;
    xlabel('Tissues'); ylabel('Sequences');
    title('Height of H3K27ac')

    % subplot(2,2,3);imagesc(permute(pcPWMp(100, :, :), [2,3,1])); colorbar;
    % subplot(2,2,4);imagesc(permute(pcPWMp(1, :, :), [2,3,1])); colorbar;


    % subplot(2,2,2);imagesc(M2(:, ind)); colorbar;
    % subplot(2,2,4);imagesc(lengths); colorbar;
end

% return the nth max of 3d matrix M
% M  - a x b x c
% maxs - a x b
function maxs = nMax(M, n)
    sortedM = sort(M, 3, 'descend');
    maxs = sortedM(:, :, n);
end

% res = 1 x k
function res = getAbsDiffMean(maxPeaks, overlaps, tissueTypeInd)
    res = mean(maxPeaks(overlaps(:, tissueTypeInd) > 0, :), 1) - mean(maxPeaks(overlaps(:, tissueTypeInd) == 0, :), 1);
    res = res - mean(maxPeaks(overlaps(:, tissueTypeInd) == 0, :), 1);
    res = abs(res);
end
function [overlaps, Xs, maxPeaks] = genData(mergedPeaksMin)
    % N = size(mergedPeaksMin.overlaps, 1);
    % N = 50;
    L = size(mergedPeaksMin.seqs, 2);
    % Xs - N x r
    % [overlaps, seqInd] = sortrows(mergedPeaksMin.overlaps(1:N,:));
    % overlaps = mergedPeaksMin.overlaps(:, [4, 12]);

    overlaps = mergedPeaksMin.overlaps(:, [1, 2]);
    % overlaps = overlaps==13;% | overlaps==14;
    size(overlaps)
    mask = mergedPeaksMin.lengths >= L;
    size(mask)
    mask = mask & (sum(overlaps > 0, 2) == 1);
    overlaps = overlaps(mask, :);
    Xs = mergedPeaksMin.seqs(mask, :);

    [overlaps, seqInd] = sortrows(overlaps);
    Xs = Xs(seqInd, :);
    size(Xs)
    size(overlaps)
    % Xs - N x L
    % Xs = mergedPeaksMin.seqs(seqInd(1:N),:);

    Xs = cat(2, Xs, fliplr(5-Xs));
    r = size(overlaps, 2); %19

    % mask = rand(size(overlaps, 1), 1)>0.50;

    % add background sequences
    % fprintf('Loading non-enhancers\n')
    % L = size(mergedPeaksMin.seqs, 2);
    % bgSeqs = matUtils.readSeq(fullfile('data', 'NEnhancers.train.seq'), L);
    % % bgSeqs = cat(2, bgSeqs, fliplr(5-bgSeqs));
    % Xs = cat(1, Xs, bgSeqs);
    % N = size(overlaps, 1);
    % overlaps = cat(1, overlaps, zeros(size(bgSeqs, 1), size(overlaps, 2)));
    % overlaps = cat(2, overlaps, [zeros(N, 1); ones(size(bgSeqs, 1), 1)]);

    fprintf('Calculating PWMs LogLikelihood\n')
    size(Xs)
    % pcPWMp = BaumWelchPWM.preComputePWMp(Xs);
    % maxPeaks = max(pcPWMp, [], 3);
    pcPWMp = BaumWelchPWM.preComputePWMpMax(Xs);
    maxPeaks = pcPWMp;


end

% function [overlaps, Xs] = genData(mergedPeaksMin)
%     mask = sum(mergedPeaksMin.overlaps, 2) == 1;
%     % Xs - N x r
%     [overlaps, seqInd] = sortrows(mergedPeaksMin.overlaps(mask,:));
%     % Xs - N x L
%     Xs = mergedPeaksMin.seqs(seqInd,:);
% end


function [mask] = removeLowHeight(overlaps, r, topPercent)
    fprintf('height\n');
    % remove low peaks
    mask = false(size(overlaps, 1), 1);
    for i = 1:r
        [~, ind] = sort(overlaps(:,i), 'descend');
        mask(ind(1:round(sum(overlaps(:,i) > 0)*topPercent))) = true;
    end
end