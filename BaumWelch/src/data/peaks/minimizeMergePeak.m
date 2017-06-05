% load('data/peaks/mergedPeaks.mat');
% [seqs, overlaps] = minimizeMergePeak(mergedPeaks);
function minimizeMergePeak(mergedPeaks)
    TOP_PEAKS_HEIGHT_PERCENT = 0.05;
    r = length(mergedPeaks(1).overlap); %19
    L = 1000;


    mergedPeaks = removeNonDistal(mergedPeaks);
    overlaps = extractOverlaps(mergedPeaks, r);
    [mergedPeaks, overlaps] = removeLowHeight(mergedPeaks, overlaps, r, TOP_PEAKS_HEIGHT_PERCENT);
    seqs = trimSeqs(mergedPeaks, L);
    [seqs, overlaps] = removeNonLetters(seqs, overlaps);

    fprintf('save\n');
    save('mergedPeaksMinimized.mat', 'seqs', 'overlaps')
end


function [seqs, overlaps] = removeNonLetters(seqs, overlaps)
    fprintf('non letters\n');
    mask = max(seqs, [], 2) <= 4;
    seqs = seqs(mask, :);
    overlaps = overlaps(mask, :);
end

function seqs = trimSeqs(mergedPeaks, L)

    fprintf('sequence\n');
    seqsCells = {mergedPeaks.seq};
    seqs = zeros(length(mergedPeaks), L);
    for i = 1:length(seqsCells)
        seq = seqsCells{i};
        center = round(length(seq)/2);
        seqs(i, :) = nt2int(seq(center-L/2+1:center+L/2));
    end
end

function [mergedPeaks, overlaps] = removeLowHeight(mergedPeaks, overlaps, r, topPercent)
    fprintf('height\n');
    % remove low peaks
    mask = false(length(mergedPeaks), 1);
    for i = 1:r
        [~, ind] = sort(overlaps(:,i), 'descend');

        mask(ind(1:round(sum(overlaps(:,i) > 0)*topPercent))) = true;
    end
    overlaps = overlaps(mask, :);
    mergedPeaks = mergedPeaks(mask);
end

function overlaps = extractOverlaps(mergedPeaks, r)

    fprintf('overlaps\n');
    overlapsFlat = [mergedPeaks.overlap];
    % N x r
    overlaps = reshape(overlapsFlat, [r, length(mergedPeaks)])';
end

function mergedPeaks = removeNonDistal(mergedPeaks)
    fprintf('Remove non distal\n');
    mergedPeaks = mergedPeaks(strcmp({mergedPeaks.class}, 'Distal'));
end