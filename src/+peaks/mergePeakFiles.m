% merges Roadmap mats to one mat file that has struct array without
% overlapping sequences, since they were all joint together and the
% overlaps vectors became from one hot to a heat map of the height of
% the peak in each tissue

% mergedPeaks fields: ['seqTo', 'peakTo', 'peakFrom', 'overlap', 'height', 'peakPos']
% withBackground - sees the background as a tissue, and takes sequences from it
% withSeq - saves a sequences actual data in file, instead only the metadata of the sequences
function [mergedPeaks, tissueEIDs, backgroundInd, genesInd] = mergePeakFiles(withBackground, withGenes, withSeq, L)
    mergedFilePath = sprintf('../data/peaks/mergedPeaks_L%db%dg%dws%d.mat', L, withBackground, withGenes, withSeq);
    fprintf('Looking for %s ...\n', mergedFilePath);
    if isfile(mergedFilePath)
        fprintf('Found %s . loading...\n', mergedFilePath);
        load(mergedFilePath);
        fprintf('Done.\n');
        return;
    end
    fprintf('Does not exist, calculating...\n');

    matFiles = peaks.beds2mats(L);
    % matFiles = peaks.beds2matsChromHMM(L);
    fprintf('Reading mat files...\n');
    [unmergedPeaks, tissueEIDs, backgroundInd, genesInd] = readMatFiles(withBackground, withGenes, matFiles);
    fprintf('Merging...\n');
    mergedPeaks = mergePeaks(unmergedPeaks, withSeq);
    save('-v7.3', mergedFilePath, 'mergedPeaks', 'tissueEIDs', 'backgroundInd', 'genesInd');
    fprintf('Saved peaks in %s\n', mergedFilePath);
end


function [unmergedPeaks, tissueEIDs, backgroundInd, genesInd] = readMatFiles(withBackground, withGenes, matFiles)
    unmergedPeaks = [];
    tissueEIDs = {};
    backgroundInd = 0;
    genesInd = 0;
    for i = 1:length(matFiles)
        filename = matFiles(i).name;
        matFilepath = fullfile(matFiles(i).folder, matFiles(i).name);
        peaks = load(matFilepath);
        fprintf('loaded mat peak data from %s\n', matFilepath);
        filenameParts = strsplit(matFiles(i).name, '_');
        tissueEID = filenameParts{1};
        if strcmp(tissueEID , 'background')
            if ~withBackground
                fprintf('skipping background\n');
                continue
            end
            backgroundInd = i;
        end
        if strcmp(tissueEID , 'genes')
            if ~withGenes
                fprintf('skipping genes\n');
                continue
            end
            genesInd = i;
        end
        if length(peaks.S) > 0
            tissueEIDs{find(peaks.S{1}.overlap)} = tissueEID;
            unmergedPeaks = [unmergedPeaks, [peaks.S{:}]];
        end
    end
end


function mergedPeaks = mergePeaks(unmergedPeaks, withSeq)

    mergedPeaks = unmergedPeaks;
    j = 1;
    for chrName = unique({unmergedPeaks.chr})
        chrMask = strcmp({unmergedPeaks.chr}, chrName{1});
        chrPeaks = unmergedPeaks(chrMask);
        [~, ind] = sort([chrPeaks.seqFrom]);
        mergedPeaks(j) = chrPeaks(ind(1));
        i = 2;
        while i <= length(ind)
            newPeak = chrPeaks(ind(i));
            oldPeak = mergedPeaks(j);
            if oldPeak.seqTo > newPeak.seqFrom
                fprintf('.')
                % merge
                if withSeq
                    oldPeak.seq = [oldPeak.seq, newPeak.seq(end - (newPeak.seqTo - oldPeak.seqTo) + 1:end)];
                    % oldPeak.seqH3K27ac = [oldPeak.seqH3K27ac, newPeak.seqH3K27ac(end - (newPeak.seqTo - oldPeak.seqTo) + 1:end)];
                    % oldPeak.seqDNase = [oldPeak.seqDNase, newPeak.seqDNase(end - (newPeak.seqTo - oldPeak.seqTo) + 1:end)];
                end
                oldPeak.seqTo = newPeak.seqTo;
                oldPeak.peakTo = max(newPeak.peakTo, oldPeak.peakTo);
                oldPeak.peakFrom = min(newPeak.peakFrom, oldPeak.peakFrom);
                % oldPeak.pos = round((oldPeak.seqFrom + newPeak.seqTo)/2) ;
                oldPeak.overlap = max(oldPeak.overlap, newPeak.overlap);
                oldPeak.height = max(oldPeak.height, newPeak.height);
                % builds an average of the peak location
                oldPeak.peakPos = ((sum(oldPeak.overlap > 0, 2) - 1) * oldPeak.peakPos + newPeak.peakPos) / sum(oldPeak.overlap > 0, 2);
                % oldPeak.min = min(oldPeak.min, newPeak.min);

                mergedPeaks(j) = oldPeak;

            else
                fprintf('O')
                j = j + 1;
                mergedPeaks(j) = newPeak;
            end
            if mod(i,100)==0
                fprintf(' %s, %d / %d\n', chrName{1}, i, length(ind));
            end
            i = i + 1;
        end
    end
    mergedPeaks = mergedPeaks(1:j);
end
