% fields of peaks.mat - ['chr', 'peakFrom', 'peakTo', 'seqFrom',
%                        'seqTo', 'peakLength', 'height', 'peakPos', 'overlap']
function typesOfCells = beds2mats()
    dbstop if error
    BEDS_DIR = '../data/peaks/processed';
    HG19_MM_DIR = '../data/peaks/mm';
    % save in dict opened hg19 fasta files as memory mapped files
    bedFiles = dir([BEDS_DIR, '/*.cleaned.narrowPeak']);
    typesOfCells = length(bedFiles);
    dict = peaks.fasta2mem();
    assert(length(dict.keys()) > 0);
    assert(length(bedFiles) > 0);
    if not(isdir('../data/peaks/mat/'))
        mkdir('../data/peaks/mat/');
    end
    for index = 1:typesOfCells
        if not(bedFiles(index).isdir)
            bedFilePath = fullfile(BEDS_DIR, bedFiles(index).name);
            nameParts = strsplit(bedFiles(index).name, '.');
            nameParts = strsplit(nameParts{1}, '-');
            name = nameParts{1};
            matPath = sprintf('../data/peaks/mat/%s.peaks.mat', name);
            if isfile(matPath)
                fprintf('file already exists, skipping. [%s]\n', matPath);
                skip;
            end
            fprintf('Converting %d / %d: [%s] %s -> %s\n', index, ...
                    typesOfCells, name, bedFilePath, matPath);
            bed2mat(index, name, bedFilePath, matPath, typesOfCells, dict);
            assert(isfile(matPath));
        end
    end
    fclose('all');
end


% cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src
function bed2mat(index, name, bedFilePath, matPath, typesOfCells, dict)
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
    if N == 0
        fprintf('no peaks found, not saving a .MAT file');
        return
    end
    % read sequences from HG19 fasta files
    fprintf(['Generating .MAT file ', name, ' (', num2str(N), ')\n']);
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
        S{newSeqId}.overlap = zeros(1, typesOfCells);
        S{newSeqId}.overlap(index) = max(heights(i), 1);
        chrLength = length(dict(S{newSeqId}.chr).Data);
        % [S{newSeqId}.seqTo, S{newSeqId}.seqFrom] = fitToL(S{newSeqId}.peakPos, L, chrLength);
        S{newSeqId}.seq = dict(S{newSeqId}.chr).Data(S{newSeqId}.seqFrom:S{newSeqId}.seqTo)';
        if mod(i, 1000) == 0
            fprintf('%%%.2f\r', 100*i/N);
        end
    end
    fprintf('\n');
    save(matPath, 'S', '-v7.3');
    fprintf('Saved bed in mat format in %s \n', matPath);
end


function dict = makeMMDict(HG19_MM_DIR)

    % save in dict opened hg19 fasta files as memory mapped files
    fprintf('Making mm dict\n');
    dict = containers.Map;
    mmFiles = dir([HG19_MM_DIR, '/*.mm']);
    assert(length(mmFiles) > 0)
    for i = 1:length(mmFiles)
        if not(mmFiles(i).isdir)
            mmFilePath = fullfile(HG19_MM_DIR, mmFiles(i).name);
            [~,filename,~] = fileparts(fullfile(HG19_MM_DIR, mmFiles(i).name));
            dict(filename) = memmapfile(mmFilePath, 'format', 'uint8');
        end
    end
end