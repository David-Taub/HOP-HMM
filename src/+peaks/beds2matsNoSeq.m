% download_and_process_all.sh

% peaks.beds2matsNoSeq()
% peaks.mergePeakFiles()
% mergedPeaks = load('../data/peaks/mergedPeaks.mat', 'mergedPeaks');
% superEnhancers = peaks.superEnhancerCaller(mergedPeaks, 10000);

% creates mats each has a cell array of only the his sequence, and have overlap one hot map that is on only in it's position
function beds2matsNoSeq()
    BEDS_DIR = '../data/peaks/processed';
    MAT_OUT_DIR = '../data/peaks/mat/';
    mkdir(MAT_OUT_DIR)
    % save in dict opened hg19 fasta files as memory mapped files
    bedFiles = dir([BEDS_DIR, '/*.cleaned.narrowPeak']);
    typesOfCells = length(bedFiles)
    for index = 1:typesOfCells
        if not(bedFiles(index).isdir)
            bedPath = fullfile(BEDS_DIR, bedFiles(index).name);
            nameParts = strsplit(bedFiles(index).name, '-');
            name = nameParts{1};
            bed2mat(index, name, bedPath, typesOfCells, MAT_OUT_DIR);
        end
    end
end

% cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src
function bed2mat(index, name, bedFilePath, typesOfCells, outDir)
    % bedPath = 'data/peaks/raw/roadmap/BrainMFLVsLiver/brain_mid_frontal_lobe/compressed/GSM773015_BI.Brain_Mid_Frontal_Lobe.H3K27ac.149.cleaned.bed';
    fprintf('Loading bed\n');
    fid = fopen(bedFilePath);
    if strcmp(name, 'genes')
        bedData = textscan(fid, '%s%d%d%*s%d%*s%d%d%d%s%s%s', 'delimiter','\t');
        chrs = bedData{1};
        peakFroms = bedData{2};
        peakTos = bedData{3};
        heights = bedData{4};
        maxPos = bedData{7};
    else
        bedData = textscan(fid,' %s%d%d%*s%d%*s%f%f%f%d', 'delimiter','\t');
        chrs = bedData{1};
        peakFroms = bedData{2};
        peakTos = bedData{3};
        heights = bedData{4};
        maxPos = bedData{8};
    end

    fclose(fid);

    N = length(chrs);
    % read sequences from HG19 fasta files
    fprintf(['Generating mat file ',name,'\n']);
    S = {};
    for i = 1:N
        newSeqId = length(S) + 1;
        S{newSeqId}.chr = chrs{i};
        S{newSeqId}.peakFrom = peakFroms(i);
        S{newSeqId}.peakTo = peakTos(i);
        % S{newSeqId}.seqFrom = peakFroms(i);
        % S{newSeqId}.seqTo = peakTos(i);
        S{newSeqId}.peakLength = S{newSeqId}.peakTo - S{newSeqId}.peakFrom;
        S{newSeqId}.height = heights(i);
        S{newSeqId}.peakPos = peakFroms(i)+maxPos(i);
        S{newSeqId}.overlap = zeros(1, typesOfCells);
        S{newSeqId}.overlap(index) = max(heights(i), 1);
        % chrLength = length(dict(S{newSeqId}.chr).Data);
        % [S{newSeqId}.seqTo, S{newSeqId}.seqFrom] = fitToL(S{newSeqId}.peakPos, L, chrLength);
        % S{newSeqId}.seq = dict(S{newSeqId}.chr).Data(S{newSeqId}.seqFrom:S{newSeqId}.seqTo)';
        if mod(i, 1000) == 0
            fprintf('%.2f\r%%', 100*i/N);
        end
    end
    fprintf('\n');

    % seqs should have 473980 sequences
    matPath = [outDir, name, '-H3k27ac.peaks.mat'];
    fprintf(['Saving mat file ', matPath, '\n']);
    save(matPath, 'S');
end