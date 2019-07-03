% merges Roadmap mats to one mat file that has struct array without
% overlapping sequences, since they were all joint together and the
% overlaps vectors became from one hot to a heat map of the height of
% the peak in each tissue

% mergedPeaks fields: ['seqTo', 'peakTo', 'peakFrom', 'overlap', 'height', 'peakPos']
% withBackground - sees the background as a tissue, and takes sequences from it
% withSeq - saves a sequences actual data in file, instead only the metadata of the sequences
function [mergedPeaks, tissueNames, backgroundInd, genesInd] = mergePeakFiles(withBackground, withGenes, withSeq, L)
    mergedFilePath = sprintf('../data/peaks/mergedPeaks_L%db%dg%dws%d.mat', L, withBackground, withGenes, withSeq);

    fprintf('Looking for %s ...\n', mergedFilePath);
    if isfile(mergedFilePath)
        fprintf('Found %s . loading...\n', mergedFilePath);
        load(mergedFilePath);
        fprintf('Done.\n');
        return;
    end
    fprintf('Does not exist, calculating...\n');
    ROADMAP_NAMES_CSV_PATH = '../data/peaks/help/full_tissue_names.csv';
    assert(isfile(ROADMAP_NAMES_CSV_PATH))

    if withSeq
        inputMatDirPath = '../data/peaks/mat';
    else
        inputMatDirPath = '../data/peaks/mat_no_seq';
    end
    matFilesCount = length(dir(inputMatDirPath)) - 2;
    if matFilesCount < 40
        fprintf('Found only %d mat files at %s. Generating mat files...\n', matFilepath, inputMatDirPath);
        peaks.beds2mats(L);
    end
    namesDict = roadmapNamesDict(ROADMAP_NAMES_CSV_PATH);
    fprintf('Reading mat files...\n');
    [unmergedPeaks, tissueNames, backgroundInd, genesInd] = readMatFiles(inputMatDirPath, withBackground, withGenes);
    tissueNames = convertNames(tissueNames, namesDict);
    fprintf('Merging...\n');
    mergedPeaks = mergePeaks(unmergedPeaks, withSeq);

    save('-v7.3', mergedFilePath, 'mergedPeaks', 'tissueNames');
    fprintf('Saved peaks in %s\n', mergedFilePath);
end

function namesDict = roadmapNamesDict(namesCSVPath)
    fid = fopen(namesCSVPath);
    csvData = textscan(fid, '%s%s', 'delimiter',',');
    fclose(fid);
    namesDict = containers.Map(csvData{1}, csvData{2});
end

function tissueNames = convertNames(tissueNames, namesDict)
    for i = 1:length(tissueNames)
        if any(strcmp(tissueNames{i}, namesDict.keys))
            key = tissueNames{i};
            tissueNames{i} = namesDict(key);
            fprintf('Tissue name found: %s -> %s\n', key, tissueNames{i});
        end
    end
end

function [unmergedPeaks, tissueNames, backgroundInd, genesInd] = readMatFiles(matDirPath, withBackground, withGenes)
    unmergedPeaks = [];
    tissueNames = {};
    backgroundInd = 0;
    genesInd = 0;
    peakFiles = dir(fullfile(matDirPath, '*.peaks.mat'));
    assert(length(peakFiles) > 0);
    for i = 1:length(peakFiles)
        filename = peakFiles(i).name;
        matFilepath = fullfile(matDirPath, filename);
        peaks = load(matFilepath);
        fprintf('loaded mat peak data from %s\n', matFilepath);
        filenameParts = strsplit(filename, '.');
        tissueName = filenameParts{1};
        if strcmp(tissueName , 'background')
            if ~withBackground
                fprintf('skipping background\n');
                continue
            end
            backgroundInd = i;
        end
        if strcmp(tissueName , 'genes')
            if ~withGenes
                fprintf('skipping genes\n');
                continue
            end
            genesInd = i;
        end
        if length(peaks.S) > 0
            tissueNames{find(peaks.S{1}.overlap)} = tissueName;
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
                    oldPeak.seq = [oldPeak.seq, newPeak.seq(end-(newPeak.seqTo-oldPeak.seqTo) + 1:end)];
                end
                oldPeak.seqTo = newPeak.seqTo;
                oldPeak.peakTo = max(newPeak.peakTo, oldPeak.peakTo);
                oldPeak.peakFrom = min(newPeak.peakFrom, oldPeak.peakFrom);
                % oldPeak.pos = round((oldPeak.seqFrom + newPeak.seqTo)/2) ;
                oldPeak.overlap = max(oldPeak.overlap, newPeak.overlap);
                oldPeak.height = max(oldPeak.height, newPeak.height);
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
