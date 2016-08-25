function [overlaps, from, to] = reader()
    suffix = '-H3K27ac.peaks.mat';
    homePath = '/cs/cbio/tommy/Enhancers/Data/GSE29184_Bing_Ren';
    fileList = dir(homePath);
    filtered = regexp({fileList.name} ,'(.+.H3K27ac.peaks.mat$)','match');
    filtered = [filtered{:}];
    cellLines = strrep(filtered, suffix, '');
    N = length(cellLines);
    chrs = {'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrM', 'chrX', 'chrY'};

    chromoLengthUB = 10^9;

    overlaps = [];
    from = [];
    to = [];
    for i = 1 : N
        i
        f = fullfile(homePath, strcat(cellLines{i}, suffix));
        load(f, 'S');
        filtered = filterPeaks(S);

        M = length(filtered);
        
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
    end
    [overlaps, from, to] = merge(overlaps, from, to);
end

function seqs = getSeqs(from, to, chrs)
    load('genome_mm9.mat')
    chrInd = floor(from / 10^9);
    chrNames = chrs(chrInd);

end

function filtered = filterPeaks(S);
    S2=[S{:}];
    heightTreshold = 0.5;
    % singular = sum(cat(1,S2.overlaps), 2) == 1;
    distal = ismember({S2.class},{'Distal'}).';
    high = ([S2.height] > S{ceil(length(S) * heightTreshold)}.height).';
    % mask = singular & distal & high;
    mask = distal & high;
    filtered = [S{mask}];
end

function [overlaps, from, to] = merge(overlaps, from, to)
    N = 19;
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

        % merge
        overlaps(1:end-1, :) = overlaps(1:end-1, :) + bsxfun(@times, overlaps(2:end, :),toMerge(2:end));
        from([toMerge(2:end);false]) = min([from([toMerge(2:end); false]), from(toMerge)], [], 2);
        to([toMerge(2:end);false]) = max([to([toMerge(2:end); false]), to(toMerge)], [], 2);

        % remove merged
        overlaps(toMerge, :) = [];
        from(toMerge) = [];
        to(toMerge) = [];
    end
end