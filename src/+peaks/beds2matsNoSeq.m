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
function bed2mat(index, name, bedPath, typesOfCells, outDir)
    % bedPath = 'data/peaks/raw/roadmap/BrainMFLVsLiver/brain_mid_frontal_lobe/compressed/GSM773015_BI.Brain_Mid_Frontal_Lobe.H3K27ac.149.cleaned.bed';
    fprintf('Loading bed\n');
    fid = fopen(bedPath);
    data = textscan(fopen(bedPath), '%s%d%d%*s%d%*s%f%f%f%d', 'delimiter','\t');
    fclose(fid);

    N = length(data{1});
    % read sequences from HG19 fasta files
    fprintf(['Generating mat file ',name,'\n']);
    S = cell(N, 1);
    for i = 1:N
        S{i}.chr = data{1}{i};
        S{i}.peakFrom = data{2}(i);
        S{i}.peakTo = data{3}(i);
        S{i}.seqFrom = data{2}(i);
        S{i}.seqTo = data{3}(i);
        S{i}.peakLength = S{i}.peakTo - S{i}.peakFrom;
        S{i}.height = data{4}(i);
        S{i}.peakPos = data{2}(i)+data{8}(i);
        S{i}.overlap = zeros(1, typesOfCells);
        S{i}.overlap(index) = 1;
        fprintf('%.2f\r%%', 100*i/N);
    end
    fprintf('\n');

    % seqs should have 473980 sequences
    matPath = [outDir, name, '-H3k27ac.peaks.mat'];
    fprintf(['Saving mat file ', matPath, '\n']);
    save(matPath, 'S');
end