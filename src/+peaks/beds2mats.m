% fields of peaks.mat - ['chr', 'peakFrom', 'peakTo', 'seqFrom',
%                        'seqTo', 'peakLength', 'height', 'peakPos', 'overlap']
function matFiles = beds2mats(L)
    dbstop if error
    BEDS_DIR = '../data/peaks/processed';
    HG19_MM_DIR = '../data/peaks/mm';
    MAX_SEQS_IN_MAT = 30000;
    % save in dict opened hg19 fasta files as memory mapped files
    bedFiles = dir([BEDS_DIR, '/*.cleaned.narrowPeak']);
    tissuesCount = length(bedFiles);
    dict = peaks.fasta2mem();
    assert(length(dict.keys()) > 0);
    assert(length(bedFiles) > 0);
    if not(isdir('../data/peaks/mat/'))
        mkdir('../data/peaks/mat/');
    end
    for index = 1:tissuesCount
        if not(bedFiles(index).isdir)
            bedFilePath = fullfile(BEDS_DIR, bedFiles(index).name);
            EIDsParts = strsplit(bedFiles(index).name, '.');
            EIDsParts = strsplit(EIDsParts{1}, '-');
            EID = EIDsParts{1};
            matPath = sprintf('../data/peaks/mat/%s_L%d_m%d.peaks.mat', EID, L, MAX_SEQS_IN_MAT);
            if isfile(matPath)
                fprintf('file already exists, skipping. [%s]\n', matPath);
                continue;
            end
            fprintf('Converting %d / %d: [%s] %s -> %s\n', index, ...
                    tissuesCount, EID, bedFilePath, matPath);
            bed2mat(index, EID, bedFilePath, matPath, tissuesCount, dict, L, MAX_SEQS_IN_MAT);
            assert(isfile(matPath));
        end
    end
    fclose('all');
    matFiles = dir(sprintf('../data/peaks/mat/*_L%d_m%d.peaks.mat', L, MAX_SEQS_IN_MAT));
    assert(length(matFiles) == length(bedFiles));
end


% cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src
function bed2mat(index, EID, bedFilePath, matPath, tissuesCount, dict, L, maxSeqsInMat)
    fprintf('Loading bed\n');
    fid = fopen(bedFilePath);
    if strcmp(EID, 'genes')
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

    % H3K27acBedGraphFilePath = sprintf("../data/peaks/processed_bedgraphs/%s-H3K27ac.enh.bedgraph", EID);
    % DNaseBedGraphFilePath = sprintf("../data/peaks/processed_bedgraphs/%s-DNase.enh.bedgraph", EID);
    % assert(isfile(H3K27acBedGraphFilePath));
    % assert(isfile(DNaseBedGraphFilePath));
    % H3K27acBedGraph = readBedGraph(H3K27acBedGraphFilePath);
    % DNaseBedGraph = readBedGraph(DNaseBedGraphFilePath);

    N = min(length(chrs), maxSeqsInMat);
    if N == 0
        fprintf('no peaks found, not saving a .MAT file');
        return
    end
    % read sequences from HG19 fasta files
    fprintf(['Generating .MAT file ', EID, ' (', num2str(N), ')\n']);
    S = {};
    for i = 1:N
        if ~any(strcmp(dict.keys(), chrs{i}))
            continue;
        end
        newSeqId = length(S) + 1;
        S{newSeqId}.chr = chrs{i};
        S{newSeqId}.peakFrom = peakFroms(i);
        S{newSeqId}.peakTo = peakTos(i);
        S{newSeqId}.seqFrom = peakFroms(i);
        S{newSeqId}.seqTo = peakTos(i);
        S{newSeqId}.peakLength = S{newSeqId}.peakTo - S{newSeqId}.peakFrom;
        S{newSeqId}.height = heights(i);
        S{newSeqId}.peakPos = peakFroms(i) + maxPos(i);
        S{newSeqId}.overlap = zeros(1, tissuesCount);
        S{newSeqId}.overlap(index) = max(heights(i), 1);
        chrLength = length(dict(S{newSeqId}.chr).Data);
        % takes sequence of length L, even if the peak is not the center is around the peak maximum
        [S{newSeqId}.seqTo, S{newSeqId}.seqFrom] = fitToL(S{newSeqId}.peakPos, L, chrLength);
        S{newSeqId}.seq = dict(S{newSeqId}.chr).Data(S{newSeqId}.seqFrom:S{newSeqId}.seqTo)';
        % S{newSeqId}.seqH3K27ac = getTrack(H3K27acBedGraph, chrs{i}, S{newSeqId}.seqFrom, S{newSeqId}.seqTo)';
        % S{newSeqId}.seqDNase = getTrack(DNaseBedGraph, chrs{i}, S{newSeqId}.seqFrom, S{newSeqId}.seqTo)';
        % S{newSeqId}.seqDnase = dict(S{newSeqId}.chr).Data(S{newSeqId}.seqFrom:S{newSeqId}.seqTo)';
        if mod(i, 1000) == 0
            fprintf('%%%.2f\r', 100*i/N);
        end
    end
    fprintf('\n');
    save(matPath, 'S', '-v7.3');
    fprintf('Saved bed in mat format in %s \n', matPath);
end


function [newTo, newFrom] = fitToL(peakPos, L, chrLength)
    center = peakPos;
    newTo = min(center + floor(L / 2), chrLength);
    newFrom = newTo - L + 1;
    if newFrom < 1
        newFrom = 1;
        newTo = L;
    end
end


function dict = makeMMDict(HG19Dir)
    % save in dict opened hg19 fasta files as memory mapped files
    fprintf('Making mm dict\n');
    dict = containers.Map;
    mmFiles = dir([HG19Dir, '/*.mm']);
    assert(length(mmFiles) > 0)
    for i = 1:length(mmFiles)
        if not(mmFiles(i).isdir)
            mmFilePath = fullfile(HG19Dir, mmFiles(i).name);
            [~,filename,~] = fileparts(fullfile(HG19Dir, mmFiles(i).name));
            dict(filename) = memmapfile(mmFilePath, 'format', 'uint8');
        end
    end
end

