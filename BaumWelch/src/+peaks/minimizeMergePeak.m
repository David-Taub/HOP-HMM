% load('data/peaks/merged/mergedPeaks.mat');
% peaks.minimizeMergePeak(mergedPeaks, 500);
function minimizeMergePeak(mergedPeaks, L)
    TOP_PEAKS_HEIGHT_PERCENT = 0.25;
    r = length(mergedPeaks(1).overlap); %19


    % mergedPeaks = removeNonDistal(mergedPeaks);
    [overlaps, lengths] = extractOverlaps(mergedPeaks, r);
    size(overlaps)
    [mergedPeaks, overlaps, lengths] = removeLowHeight(mergedPeaks, overlaps, lengths, r, TOP_PEAKS_HEIGHT_PERCENT);
    size(overlaps)
    seqs = trimSeqs(mergedPeaks, L);
    [seqs, overlaps, lengths] = removeNonLetters(seqs, overlaps, lengths);
    size(overlaps)

    fprintf('save\n');
    save('data/peaks/roadmap/mergedPeaksMinimized.mat', 'seqs', 'overlaps', 'lengths')
end


function [seqs, overlaps, lengths] = removeNonLetters(seqs, overlaps, lengths)
    fprintf('non letters\n');
    mask = max(seqs, [], 2) <= 4;
    seqs = seqs(mask, :);
    overlaps = overlaps(mask, :);
    lengths = lengths(mask);
end

function seqs = trimSeqs(mergedPeaks, L)

    fprintf('sequence\n');
    seqsCells = {mergedPeaks.seq};
    seqs = zeros(length(mergedPeaks), L);
    for i = 1:length(seqsCells)
        seq = seqsCells{i};
        center = round(length(seq)/2);
        % seqs(i, :) = nt2int(seq(center-L/2+1:center+L/2));
        start_pos = center- L / 2 +1;
        end_pos = start_pos + L;
        seqs(i, :) = seq(start_pos:end_pos);
    end
end

function [mergedPeaks, overlaps, lengths] = removeLowHeight(mergedPeaks, overlaps, lengths, r, topPercent)
    fprintf('height\n');
    % remove low peaks
    mask = false(length(mergedPeaks), 1);
    for i = 1:r
        [~, ind] = sort(overlaps(:,i), 'descend');

        mask(ind(1:round(sum(overlaps(:,i) > 0)*topPercent))) = true;
    end
    overlaps = overlaps(mask, :);
    mergedPeaks = mergedPeaks(mask);
    lengths = lengths(mask);
end

function [overlaps, lengths] = extractOverlaps(mergedPeaks, r)

    fprintf('overlaps\n');
    overlapsFlat = [mergedPeaks.overlap];
    lengths = [mergedPeaks.peakLength]';
    % N x r
    overlaps = reshape(overlapsFlat, [r, length(mergedPeaks)])';
end

function mergedPeaks = removeNonDistal(mergedPeaks)
    fprintf('Remove non distal\n');
    mergedPeaks = mergedPeaks(strcmp({mergedPeaks.class}, 'Distal'));
end