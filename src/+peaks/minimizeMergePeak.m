% fields of minimizeMergePeak - ['mergedPeaksMin.seqs', 'mergedPeaksMin.overlaps', 'mergedPeaksMin.peakLengths', 'mergedPeaksMin.tissueNames']
function mergedPeaksMin = minimizeMergePeak(topPercent, doEnhSpecific, withBackground, withGenes,...
                                            seqsPerTissue, L, peakMinL, peakMaxL, tissueList, minSamplesCount)
    minimizedMergedFilePath = sprintf('../data/peaks/mergedPeaksMinimized_L%db%dg%dp%des%dspt%dmin%dmax%dsc%dT%s.mat', L, ...
                                      withBackground, withGenes, floor(100 * topPercent), doEnhSpecific, seqsPerTissue, ...
                                      peakMinL, peakMaxL, minSamplesCount, sprintf('%d', tissueList));
    fprintf('Looking for %s ...\n', minimizedMergedFilePath);
    if isfile(minimizedMergedFilePath)
        fprintf('Found %s . loading...\n', minimizedMergedFilePath);
        load(minimizedMergedFilePath);
        fprintf('Done.\n');
        return
    end
    fprintf('Does not exist, calculating...\n');
    [mergedPeaks, tissueEIDs, backgroundInd, genesInd] = peaks.mergePeakFiles(withBackground, withGenes, true, L);
    mergedPeaksMin.backgroundInd = backgroundInd;
    mergedPeaksMin.genesInd = genesInd;
    mergedPeaksMin.tissueEIDs = tissueEIDs;
    mergedPeaksMin.tissueNames = misc.EIDsToTissueNames(tissueEIDs);
    mergedPeaksMin = extractFields(mergedPeaks, mergedPeaksMin);

    fprintf('Seqs %d tissues [%s]\n', size(mergedPeaksMin.overlaps, 1), sprintf('%d ', sum(mergedPeaksMin.overlaps > 0, 1)));
    % the merging of the peaks result in longer peaks (unless they are in the exact same location)
    mergedPeaksMin.seqs = extractSeqs(mergedPeaks, L);
    mergedPeaksMin.samplesCount = zeros(size(mergedPeaksMin.seqs, 1), 1);

    mergedPeaksMin = removeLowHeight(mergedPeaksMin, topPercent);
    fprintf('Seqs %d tissues [%s]\n', size(mergedPeaksMin.overlaps, 1), sprintf('%d ', sum(mergedPeaksMin.overlaps > 0, 1)));

    mergedPeaksMin = removeNonLetters(mergedPeaksMin);
    fprintf('Seqs %d tissues [%s]\n', size(mergedPeaksMin.overlaps, 1), sprintf('%d ', sum(mergedPeaksMin.overlaps > 0, 1)));

    if (peakMaxL > 0) & (peakMinL > 0)
        mergedPeaksMin = removeByLength(mergedPeaksMin, peakMinL, peakMaxL)
        fprintf('Seqs %d tissues [%s]\n', size(mergedPeaksMin.overlaps, 1), sprintf('%d ', sum(mergedPeaksMin.overlaps > 0, 1)));
    end
    if doEnhSpecific
        mergedPeaksMin = removeNonSpecific(mergedPeaksMin);
        fprintf('Seqs %d tissues [%s]\n', size(mergedPeaksMin.overlaps, 1), sprintf('%d ', sum(mergedPeaksMin.overlaps > 0, 1)));
    end
    if length(tissueList) > 1 | withBackground | withGenes
        mergedPeaksMin = removeByTissueList(mergedPeaksMin, tissueList);
        fprintf('Seqs %d tissues [%s]\n', size(mergedPeaksMin.overlaps, 1), sprintf('%d ', sum(mergedPeaksMin.overlaps > 0, 1)));
    end
    mergedPeaksMin = removeBySampleCount(mergedPeaksMin, minSamplesCount);
    fprintf('Seqs %d tissues [%s]\n', size(mergedPeaksMin.overlaps, 1), sprintf('%d ', sum(mergedPeaksMin.overlaps > 0, 1)));
    if seqsPerTissue > 0
        mergedPeaksMin = balanceOverlaps(mergedPeaksMin, seqsPerTissue);
        fprintf('Seqs %d tissues [%s]\n', size(mergedPeaksMin.overlaps, 1), sprintf('%d ', sum(mergedPeaksMin.overlaps > 0, 1)));
    end
    % outFilepath = '../data/peaks/mergedPeaksMinimized.mat';
    save(minimizedMergedFilePath, '-v7.3', 'mergedPeaksMin');
    fprintf('Saved peaks in %s\n', minimizedMergedFilePath);
end



function mergedPeaksMin = removeByTissueList(mergedPeaksMin, tissueList);
    fprintf('tissue list\n');
    if mergedPeaksMin.backgroundInd > 0
        tissueList = [tissueList, mergedPeaksMin.backgroundInd];
    end
    if mergedPeaksMin.genesInd > 0
        tissueList = [tissueList, mergedPeaksMin.genesInd];
    end
    mask = sum(mergedPeaksMin.overlaps(:, tissueList) > 0, 2) > 0;
    mergedPeaksMin = reduceData(mask, mergedPeaksMin);
    mergedPeaksMin.overlaps = mergedPeaksMin.overlaps(:, tissueList);
    mergedPeaksMin.tissueNames = mergedPeaksMin.tissueNames(tissueList);
    mergedPeaksMin.tissueEIDs = mergedPeaksMin.tissueEIDs(tissueList);
    mergedPeaksMin.tissueList = tissueList;
end


function mergedPeaksMin = removeNonSpecific(mergedPeaksMin)
    fprintf('enhancer specific\n');
    mask = sum(mergedPeaksMin.overlaps > 0, 2) == 1;
    mergedPeaksMin = reduceData(mask, mergedPeaksMin);
end


function mergedPeaksMin = removeBySampleCount(mergedPeaksMin, minSamplesCount)
    fprintf('filter by samples count\n');
    mergedPeaksMin = countSamplesNumber(mergedPeaksMin);
    mask = mergedPeaksMin.samplesCount > minSamplesCount;
    if (mergedPeaksMin.backgroundInd > 0) | (mergedPeaksMin.genesInd > 0)
        mask = mask | mergedPeaksMin.overlaps(:, end) > 0
        if (mergedPeaksMin.backgroundInd > 0) & (mergedPeaksMin.genesInd > 0)
            mask = mask | mergedPeaksMin.overlaps(:, end - 1) > 0
        end
    end
    mergedPeaksMin = reduceData(mask, mergedPeaksMin);
end


function mergedPeaksMin = removeByLength(mergedPeaksMin, peakMinL, peakMaxL)
    fprintf('peak length\n');
    mask = peakMinL < mergedPeaksMin.peakLengths < peakMaxL;
    mergedPeaksMin = reduceData(mask, mergedPeaksMin);
end


function mergedPeaksMin = balanceOverlaps(mergedPeaksMin, seqsPerTissue)
    fprintf('balancing mergedPeaksMin.seqs count per tissue\n');
    mask = false(size(mergedPeaksMin.seqs, 1), 1);
    foundMaxSeqsPerTissue = sum(mergedPeaksMin.overlaps > 0, 1);
    foundMaxSeqsPerTissue = min(foundMaxSeqsPerTissue, [], 2);
    for i = 1:size(mergedPeaksMin.overlaps, 2)
        [vals, inds] = sort(mergedPeaksMin.overlaps(:, i), 1, 'descend');
        assert(vals(foundMaxSeqsPerTissue) > 0);
        mask(inds(1:foundMaxSeqsPerTissue)) = true;
    end
    mergedPeaksMin = reduceData(mask, mergedPeaksMin);
end


function mergedPeaksMin = removeNonLetters(mergedPeaksMin)
    fprintf('non letters\n');
    mask = max(mergedPeaksMin.seqs, [], 2) <= 4;
    mergedPeaksMin = reduceData(mask, mergedPeaksMin);
end


% remove sequences with peaks values that are not inside the top percentage in any tissue
function mergedPeaksMin = removeLowHeight(mergedPeaksMin, topPercent)
    fprintf('height\n');
    numerOfTissues = size(mergedPeaksMin.overlaps, 2);
    % remove low peaks
    mask = false(length(mergedPeaksMin.seqs), 1);
    for i = 1:numerOfTissues
        [vals, ind] = sort(mergedPeaksMin.overlaps(:, i), 'descend');
        amountToKeep = round(sum(mergedPeaksMin.overlaps(:, i) > 0) * topPercent);
        if vals(1) == vals(end)
            % not an enhancer height peak, either background or gene peak
            amountToKeep = sum(mergedPeaksMin.overlaps(:, i) > 0);
        end
        mask(ind(1:amountToKeep)) = true;
    end
    mergedPeaksMin = reduceData(mask, mergedPeaksMin);
end


function seqs = extractSeqs(mergedPeaks, L)
    fprintf('extract sequence\n');
    seqsCells = {mergedPeaks.seq};
    seqs = zeros(length(mergedPeaks), L);
    for i = 1:length(seqsCells)
        seq = seqsCells{i};
        center = round(length(seq) / 2);
        % seqs(i, :) = nt2int(seq(center-L/2+1:center+L/2));
        startPeakPos = center - L / 2 + 1;
        endPeakPos = startPeakPos + L - 1 ;
        seqs(i, :) = seq(startPeakPos:endPeakPos);
    end
end


% TODO: get seq start end and chrs
function minimizeMergePeak = extractFields(mergedPeaks, minimizeMergePeak)
    fprintf('fields extraction\n');
    numerOfTissues = length(mergedPeaks(1).overlap);
    overlapsFlat = [mergedPeaks.overlap];
    minimizeMergePeak.peakLengths = [mergedPeaks.peakLength]';
    minimizeMergePeak.chrs = {mergedPeaks.chr};
    minimizeMergePeak.starts = [mergedPeaks.seqFrom]';
    % N x numerOfTissues
    minimizeMergePeak.overlaps = reshape(overlapsFlat, [numerOfTissues, length(mergedPeaks)])';
end


function mergedPeaksMin = reduceData(mask, mergedPeaksMin);
    fprintf('keeping %.2f%% of mergedPeaksMin.seqs\n', 100 * mean(mask, 1));
    mergedPeaksMin.seqs = mergedPeaksMin.seqs(mask, :, :);
    mergedPeaksMin.overlaps = mergedPeaksMin.overlaps(mask, :);
    mergedPeaksMin.peakLengths = mergedPeaksMin.peakLengths(mask);
    mergedPeaksMin.starts = mergedPeaksMin.starts(mask);
    mergedPeaksMin.chrs = mergedPeaksMin.chrs(mask);
    mergedPeaksMin.samplesCount = mergedPeaksMin.samplesCount(mask);
end


function mergedPeaksMin = countSamplesNumber(mergedPeaksMin)
    fprintf('counts samples\n');
    trackNames = {'H3K27ac', 'DNase'};
    bedGraphs = misc.readAllBedGraphs(mergedPeaksMin.tissueEIDs, trackNames);
    [N, L] = size(mergedPeaksMin.seqs);
    mergedPeaksMin.samplesCount = zeros(N, 1);
    for i = 1:N
        chr = mergedPeaksMin.chrs{i};
        from = mergedPeaksMin.starts(i);
        to = mergedPeaksMin.starts(i) + L;
        mergedPeaksMin.samplesCount(i) = countSamples(bedGraphs, chr, from, to);
        if mod(i, 100) == 0
            fprintf('count sequence %d, %s: %d-%d [%d]\n', i, chr, from, to, mergedPeaksMin.samplesCount(i));
        end
    end
end


function samplesCount = countSamples(bedGraphs, trackChr, trackFrom, trackTo)
    samplesCount = inf;
    for bedGraph = [bedGraphs{:}]
        mask = strcmp(bedGraph.chrs, trackChr) & (bedGraph.tos >= trackFrom) & (bedGraph.froms <= trackTo);
        % linear interp
        samplesCount = min(samplesCount, sum(mask));
    end
end