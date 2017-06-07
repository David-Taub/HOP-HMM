% checks correlation between the existence of p300 peak and its TF binding sites
% presence. the binding sights presence are evaluated with PWMs.


% usage:
% cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src
% mergedPeaksMin = load(fullfile('data', 'peaks', 'mergedPeaksMinimized.mat'));
% PC_PWM_PROBABILITY_FILE = fullfile('data', 'precomputation', 'pcPWMp.mat');
% delete(PC_PWM_PROBABILITY_FILE);
% mainPWMCor(mergedPeaksMin);
function [maxPeaks, overlaps] = mainPWMCor(mergedPeaksMin)
    close all;
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
    tic
    ranks = getRanks(maxPeaks, overlaps, r, k);
    toc
    fprintf('showing results 1\n')
    showData(maxPeaks, ranks, overlaps);
    fprintf('showing results 2\n')
    showBestSep(ranks, maxPeaks, overlaps, r);
    save(fullfile('data', 'outMainPWMCor.mat'));
end

% maxPeaks - N x k x M
% overlaps - N x r
% ranks - r x k x M
function [ranks] = getRanks(maxPeaks, overlaps, r, k)
    ranks = zeros(r, k);
    for i = 1:r
        % 1 x k
        for j = 1:k
            peaksIndicatorTFPos = maxPeaks(overlaps(:, i) > 0, j);
            peaksIndicatorTFNeg = maxPeaks(overlaps(:, i) == 0, j);
            if mean(peaksIndicatorTFPos) > mean(peaksIndicatorTFNeg);
                ranks(i, j) = getAucRoc(peaksIndicatorTFPos, peaksIndicatorTFNeg, false);
            else
                ranks(i, j) = getAucRoc(peaksIndicatorTFNeg, peaksIndicatorTFPos, false);
            end
            % [~, ranks(i, j)] = kstest2(peaksIndicatorTFPos, peaksIndicatorTFNeg);
            % ranks(i, j) = 1 - ranks(i, j);

            % [ranks(i, j), ~] = ranksum(peaksIndicatorTFPos, peaksIndicatorTFNeg);
            % ranks(i, j) = 1 - ranks(i, j);
            % ranks(i, j) = matUtils.emdTest(peaksIndicatorTFPos, peaksIndicatorTFNeg);
            fprintf('\r%d / %d. %d / %d', i, r, j, k);
        end
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
    tissueIndicator = zeros(r, 1);
    for i =1:r
        % 1 x k
        [bestRank, tissueIndicator(i)] = max(ranks(i, :), [], 2);
        peaksIndicatorTFPos = maxPeaks(overlaps(:, i) > 0, tissueIndicator(i));
        peaksIndicatorTFNeg = maxPeaks(overlaps(:, i) == 0, tissueIndicator(i));
        % subplot(4,5,i);
        figure;
        subplot(1,2,1);
        h = histogram(peaksIndicatorTFPos, 50, 'Normalization', 'probability');
        hold on;
        histogram(peaksIndicatorTFNeg, h.BinEdges, 'Normalization', 'probability');
        title(['PWM ', int2str(tissueIndicator(i)), ' LogLikes. Rate: ', num2str(bestRank)]);
        legend('Enhancers from tissue', 'Enhancers from all other tissues')
        subplot(1,2,2);
        if mean(peaksIndicatorTFPos) > mean(peaksIndicatorTFNeg);
            getAucRoc(peaksIndicatorTFPos, peaksIndicatorTFNeg, true);
        else
            getAucRoc(peaksIndicatorTFNeg, peaksIndicatorTFPos, true);
        end
        drawnow
        % keyboard
    end
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

function showData(maxPeaks, meansDiff, overlaps)

    [maxPeaks, motifInd] = matUtils.clustMatRows(maxPeaks');
    [maxPeaks, seqInd] = matUtils.clustMatRows(maxPeaks');
    subplot(1,3,1);imagesc(maxPeaks); colorbar;
    xlabel('motif'); ylabel('sequences');
    title('max peak of motif in sequence');


    [meansDiff, tissueInd]  = matUtils.clustMatRows(meansDiff);
    meansDiff = meansDiff(:, motifInd);
    subplot(1,3,2);imagesc(meansDiff); colorbar;
    xlabel('motif'); ylabel('tissue');
    title('abs diff between mean PWM LogLikes (pos tissue vs neg tissue)')


    overlaps = overlaps(seqInd, :);
    overlaps = overlaps(:, tissueInd);
    subplot(1,3,3);imagesc(overlaps); colorbar;
    xlabel('tissues'); ylabel('sequences');
    title('Overlaps')

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

    % Xs - N x r
    % [overlaps, seqInd] = sortrows(mergedPeaksMin.overlaps(1:N,:));
    [overlaps, seqInd] = sortrows(mergedPeaksMin.overlaps);
    overlaps = overlaps(:, [5, 12]);
    inds = sum(overlaps > 0, 2) == 1;
    overlaps = overlaps(inds, :);
    % Xs - N x L
    % Xs = mergedPeaksMin.seqs(seqInd(1:N),:);
    Xs = mergedPeaksMin.seqs(seqInd(inds),:);
    Xs = cat(2, Xs, fliplr(5-Xs));
    fprintf('calculating PWMs LogLikelihood\n')
    r = size(overlaps, 2); %19
    mask = removeLowHeight(overlaps, r, 0.15);
    overlaps = overlaps(mask, :);
    Xs = Xs(mask, :);
    % pcPWMp = BaumWelchPWM.preComputePWMp(Xs);
    % maxPeaks = sum(pcPWMp, 3);
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