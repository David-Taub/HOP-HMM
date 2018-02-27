% merges Roadmap mats to one mat file that has struct array without
% overlapping sequences, since they were all joint together and the
% overlaps vectors became from one hot to a heat map of the height of
% the peak in each tissue

% download_and_process_all.sh

% first we choose the PWMs that are most different between the following tissues,
% including background and genes:

% peaks.beds2mats(500)
% peaks.mergePeakFiles(true, true)
% load('../data/peaks/mergedPeaks.mat');
% minimizeMergePeak(mergedPeaks, L, tissueNames)
% mergedPeaksMin = load('../data/peaks/mergedPeaksMinimized.mat')
% chooseBestPWMs(mergedPeaksMin, [10, 20, 30, 45, 46])

% then we find superEnhancers, and try to learn them
% peaks.beds2matsNoSeq()
% peaks.mergePeakFiles(false, false)
% load('../data/peaks/mergedPeaksNoBackground.mat');
% superEnhancers = peaks.superEnhancerCaller(mergedPeaks, 10000, [10, 20, 30]);
% mainRealData(superEnhancers, 4, 40);

function mergedPeaks = mergePeakFiles(withBackground, with_seq)
    dbstop if error
    if with_seq
        INPUT_MAT_DIR_PATH = '../data/peaks/mat';
    else
        INPUT_MAT_DIR_PATH = '../data/peaks/mat_no_seq';
    end
    OUT_FILE_PATH = '../data/peaks/mergedPeaks.mat';
    if ~withBackground
        OUT_FILE_PATH = '../data/peaks/mergedPeaksNoBackground.mat';
    end
    ROADMAP_NAMES_CSV_PATH = '../data/peaks/help/full_tissue_names.csv';
    namesDict = roadmapNamesDict(ROADMAP_NAMES_CSV_PATH);
    fprintf('Reading mat files...\n')
    [unmergedPeaks, tissueNames] = readMatFiles(INPUT_MAT_DIR_PATH, withBackground);
    tissueNames = convertNames(tissueNames, namesDict);
    fprintf('Merging...\n')
    mergedPeaks = mergePeaks(unmergedPeaks, with_seq);
    fprintf('Saving...\n')

    save('-v7.3', OUT_FILE_PATH, 'mergedPeaks', 'tissueNames');
    fprintf('Done\n')
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
            tissueNames{i} = namesDict(tissueNames{i});
        end
    end
end

function [unmergedPeaks, tissueNames] = readMatFiles(matDirPath, withBackground)
    unmergedPeaks = [];
    tissueNames = {};
    peakFiles = dir(fullfile(matDirPath, '*.peaks.mat'));
    for i = 1:length(peakFiles)
        filename = peakFiles(i).name;
        peaks = load(fullfile(matDirPath, filename));
        filenameParts = strsplit(filename, '.');
        tissueName = filenameParts{1};
        if strcmp(tissueName , 'background') && ~withBackground
            continue
        end
        tissueNames{find(peaks.S{1}.overlap)} = tissueName;
        unmergedPeaks
        unmergedPeaks = [unmergedPeaks, [peaks.S{:}]];
    end
end

function mergedPeaks = mergePeaks(unmergedPeaks, with_seq)

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
                if with_seq
                    oldPeak.seq = [oldPeak.seq, newPeak.seq(end-(newPeak.seqTo-oldPeak.seqTo) + 1:end)];
                end
                oldPeak.seqTo = newPeak.seqTo;
                oldPeak.peakTo = max(newPeak.peakTo, oldPeak.peakTo);
                oldPeak.peakFrom = min(newPeak.peakFrom, oldPeak.peakFrom);
                % oldPeak.pos = round((oldPeak.seqFrom + newPeak.seqTo)/2) ;
                oldPeak.overlap = max(oldPeak.overlap, newPeak.overlap);
                oldPeak.height = max(oldPeak.height, newPeak.height);
                oldPeak.peakPos = mean([oldPeak.peakPos, newPeak.peakPos], 2);
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
