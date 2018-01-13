
% cd /cs/stud/boogalla/projects/CompGenetics/BaumWelch/src
% peaks.readRoadMapCSV(1000);
function readRoadMapCSV(L)
    CSV_PATH = 'data/peaks/raw/roadmap/ms.csv';
    HG19_MM_DIR = 'data/peaks/raw/roadmap/hg19/mm';

    fprintf('Loading roadmap CSV\n');
    fid = fopen(CSV_PATH);
    names = textscan(fid, repmat('%q', [1, 129]),1 , 'delimiter',',');
    data = textscan(fid, ['"%[^:]:%d-%d"', repmat('%d', [1, 127])], 'HeaderLines',1, 'delimiter',',');
    fclose(fid);

    names = names(2:end-1);
    N = length(data{1});
    overlaps = data(:,4:end);
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

    % read sequences from HG19 fasta files
    fprintf('Generating mat file\n');
    seqs = cell(N, 1);
    lengths = zeros(N, 1);
    for i = 1:N
        chr = data{1}{i};
        from = data{2}(i);
        to = data{3}(i);
        lengths(i) = to - from;
        [to, from] = fitToL(to, from, L);
        seqs{i} = dict(chr).Data(from:to);
        fprintf('%.2f\r%%', 100*i/N);
    end
    fprintf('\n');

    overlaps = [overlaps{:}];
    seqs = [seqs{:}]';

    toRemove = any(seqs>4, 2);
    seqs(toRemove, :) = [];
    overlaps(toRemove, :) = [];
    lengths(toRemove, :) = [];

    % seqs should have 473980 sequences
    fprintf('Saving mat file\n');
    save('data/RoadmapEnhancers.mat', 'seqs', 'overlaps', 'names', 'lengths');
end

function [newTo, newFrom] = fitToL(to, from, L)
    center = round((to + from) / 2);
    newTo = center + floor(L / 2);
    newFrom = center - floor(L / 2) + 1;
end
