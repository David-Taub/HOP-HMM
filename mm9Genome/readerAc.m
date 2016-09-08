function [seqs, overlaps, from, to] = readerAc()
    close all;
    suffix = '-H3K27ac.peaks.mat';

    % homePath = '/cs/cbio/tommy/Enhancers/Data/GSE29184_Bing_Ren';
    H3K4me3 = '/cs/cbio/david/data/mat/H3K4me3';
    H3K4me1 = '/cs/cbio/david/data/mat/H3K4me1';
    H3K27ac = '/cs/cbio/david/data/mat/H3K27ac';
    % load('/cs/cbio/tommy/Enhancers/Data/GSE29184_Bing_Ren/spleen-H3K27ac.peaks.mat');
    fileList = dir(H3K27ac);
    fileList = {fileList.name};
    fileList = {fileList{3:end}};
    N = length(fileList);
    chrs = {'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrM', 'chrX', 'chrY'};
    chromoLengthUB = 10^9;

    tic
    overlaps = [];
    from = [];
    to = [];
    for i = 1 : N
        % i
        fileList{i}
        load(fullfile(H3K27ac, fileList{i}), 'S');
        acPeaks = S;
        filtered = filterPeaks(acPeaks);

        M = length(filtered);
        overlapsAdd = [];
        overlapsAdd(1:M, 1:i-1) = 0;
        overlapsAdd(1:M, i) = 1;
        overlapsAdd(1:M, i+1:N) = 0;
        [~,chrInd] = ismember({filtered.chr}, chrs);
        % bit Arab, to distinguish the addresses of peaks in two different
        % chromosomes, I add the chromosome index, times 1 billion to the address
        fromAdd = [filtered.from] + (chromoLengthUB .* chrInd);
        toAdd = [filtered.to] + (chromoLengthUB .* chrInd);

        overlaps = cat(1, overlaps, overlapsAdd);
        from = cat(1, from, fromAdd.');
        to = cat(1, to, toAdd.');
        length(toAdd)
    end
    [overlaps, from, to] = merge(overlaps, from, to);
    load('/cs/cbio/tommy/Enhancers/Data/genome_mm9.mat');
    seqs = getSeqs(from, to, genome);
    length(seqs)
    % peaksPath = strcat('peaks', num2str(k), '.mat');
    % save(peaksPath, 'seqs', 'overlaps', 'from', 'to');
    save('/cs/stud/boogalla/projects/CompGenetics/BaumWelch/peaks.mat', 'seqs', 'overlaps', 'from', 'to');
    toc
end

function seqs = getSeqs(from, to, genome)
    chromoLengthUB = 10^9;
    chrs = fieldnames(genome);
    chrInd = floor(from / chromoLengthUB);
    from = mod(from, chromoLengthUB);
    to = mod(to, chromoLengthUB);
    chrNames = chrs(chrInd);
    for i = 1:length(chrNames)
        fullChromosome = getfield(genome, chrNames{i});
        seqs{i} = fullChromosome(from(i):to(i));
    end
end

function filtered = filterPeaks(acPeaks)
    filtered = acPeaks;
    acS = [filtered{:}];
    
    % distance
    minDistance = 3000;
    distal = abs([acS.gTSSdist]) > minDistance;
    filtered = {filtered{distal}};

    % height
    acS = [filtered{:}];

    [vals, ord] = sort([acS.height], 'descend');
    filtered = {filtered{ord(1:15000)}};
    % filtered = [filtered{ord(vals > 30)}];
    % filtered = {filtered{ord(1:ceil(length(acS) * 0.03  ))}};
    filtered = [filtered{:}];
end

function [overlaps, from, to] = merge(overlaps, from, to)
    [~, ord] = sort(from);
    % sort
    from = from(ord);
    to = to(ord);
    overlaps = overlaps(ord, :);
    
    M = 0;
    while M ~= size(from, 1)
        M = size(from, 1);
        % center
        centers = (from + to) ./ 2;

        % mark to merge
        toMerge = false(M, 1);
        toMerge(2:end) = from(1:end-1) < centers(2:end);
        toMerge(2:end) = toMerge(2:end) & centers(2:end) < to(1:end-1);
        % a = sum(toMerge, 1);
        toMerge(2:end) = diff(toMerge) > 0;
        % b = sum(toMerge, 1);
        % [M, a, b]

        % mergeKDF
        overlaps(1:end-1, :) = overlaps(1:end-1, :) + ...
         bsxfun(@times, overlaps(2:end, :),toMerge(2:end));
        from([toMerge(2:end);false]) = min([from([toMerge(2:end); false]), from(toMerge)], [], 2);
        to([toMerge(2:end);false]) = max([to([toMerge(2:end); false]), to(toMerge)], [], 2);

        % remove merged
        overlaps(toMerge, :) = [];
        from(toMerge) = [];
        to(toMerge) = [];
    end
end