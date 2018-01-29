% peaks.beds2mats(500)
% peaks.mergePeakFiles()
% load('../data/peaks/mergedPeaks.mat');
% peaks.minimizeMergePeak(mergedPeaks, 500);
% creates mats each has a cell array of only the his sequence, and have overlap one hot map that is on only in it's position

function beds2mats(L)
    dbstop if error
    BEDS_DIR = '../data/peaks/processed';
    % save in dict opened hg19 fasta files as memory mapped files
    bedFiles = dir([BEDS_DIR, '/*.cleaned.narrowPeak']);
    typesOfCells = length(bedFiles)
    dict = peaks.fasta2mem();
    for index = 1:typesOfCells
        if not(bedFiles(index).isdir)
            bedPath = fullfile(BEDS_DIR, bedFiles(index).name);
            nameParts = strsplit(bedFiles(index).name, '-');
            name = nameParts{1};
            bed2mat(index, name, bedPath, typesOfCells, L, dict);
        end
    end
    fclose('all')
end

% cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src
function bed2mat(index, name, bedPath, typesOfCells, L, dict)
    fprintf('Loading bed\n');
    fid = fopen(bedPath);
    bedData = textscan(fopen(bedPath), '%s%d%d%*s%d%*s%f%f%f%d', 'delimiter','\t');
    fclose(fid);

    N = length(bedData{1});
    % read sequences from HG19 fasta files
    fprintf(['Generating mat file ',name,'\n']);
    S = {};
    for i = 1:N
        if ~any(strcmp(dict.keys(), bedData{1}{i}))
            continue;
        end
        newSeqId = length(S) + 1;
        S{newSeqId}.chr = bedData{1}{i};
        S{newSeqId}.peakFrom = bedData{2}(i);
        S{newSeqId}.peakTo = bedData{3}(i);
        S{newSeqId}.seqFrom = bedData{2}(i);
        S{newSeqId}.seqTo = bedData{3}(i);
        S{newSeqId}.peakLength = S{newSeqId}.peakTo - S{newSeqId}.peakFrom;
        S{newSeqId}.height = bedData{4}(i);
        S{newSeqId}.peakPos = bedData{2}(i)+bedData{8}(i);
        S{newSeqId}.overlap = zeros(1, typesOfCells);
        S{newSeqId}.overlap(index) = 1;
        chrLength = length(dict(S{newSeqId}.chr).Data);
        [S{newSeqId}.seqTo, S{newSeqId}.seqFrom] = fitToL(S{newSeqId}.peakPos, L, chrLength);
        S{newSeqId}.seq = dict(S{newSeqId}.chr).Data(S{newSeqId}.seqFrom:S{newSeqId}.seqTo)';
        if mod(i, 1000) == 0
            fprintf('%.2f\r%%', 100*i/N);
        end
    end
    fprintf('\n');

    % seqs should have 473980 sequences

    matPath = ['../data/peaks/mat/', name, '.peaks.mat'];
    fprintf(['Saving mat file ', matPath, '\n']);
    save(matPath, 'S', '-v7.3');
end

% function [newTo, newFrom] = fitToL(to, from, L)
function [newTo, newFrom] = fitToL(peakPos, L, chrLength)
    center = peakPos;
    newTo = min(center + floor(L / 2), chrLength);
    newFrom = newTo - L + 1;
    if newFrom < 1
        newFrom = 1;
        newTo = L;
    end
end

function dict = makeMMDict()
    HG19_MM_DIR = '../data/peaks/mm';
    % save in dict opened hg19 fasta files as memory mapped files
    fprintf('Making mm dict\n');
    dict = containers.Map;
    mmFiles = dir([HG19_MM_DIR, '/*.mm']);
    for i = 1:length(mmFiles)
        if not(mmFiles(i).isdir)
            mmFilePath = fullfile(HG19_MM_DIR, mmFiles(i).name);
            [~,filename,~] = fileparts(fullfile(HG19_MM_DIR, mmFiles(i).name));
            dict(filename) = memmapfile(mmFilePath, 'format', 'uint8');
        end
    end
end