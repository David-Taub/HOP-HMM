function bgSeqs = sampleBG(N, L, genes, genome, enhancersFrom, enhancersTo);
    tic
    close all;
    chrs = fieldnames(genome);
    % chrs = {'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrM', 'chrX', 'chrY'};
    chromoLengthUB = 10^9;
    gap = 7000;
    [genesFrom, genesTo] = readGenes(genes, chromoLengthUB, chrs);
    from = [enhancersFrom; genesFrom];
    to = [enhancersTo; genesTo];
    
    % add "saftey" gap before and after enhancer \ gene
    t = mod(from, chromoLengthUB) - gap;
    t < 0 == 0;
    from = floor(from ./ chromoLengthUB) .* chromoLengthUB + t;
    to = to + gap;

    % merge overlaping gene and enhancers (TODO: exist?)
    [from, to] = merge(from, to);
    bgFrom = cat(1, to, [1 : length(chrs)].' .* chromoLengthUB + 1);
    chrsLength = cellfun(@length, struct2cell(genome));
    bgTo = cat(1, from, chrsLength + [1 : length(chrs)].' .* chromoLengthUB);

    % remove BG regions which are too small
    ind = (bgTo - bgFrom < L);
    bgFrom(ind) = [];
    bgTo(ind) = [];

    % sample from background regions
    bgSeqs = sampleSeqs(bgFrom, bgTo, genome, chrs, chromoLengthUB, N, L);
    save('bg.mat', 'bgSeqs')
    toc
end

function seqs = sampleSeqs(from, to, genome, chrs, chromoLengthUB, N, L)
    seqs = [];
    while size(seqs, 1) < N
        addSeqsAmount = N - size(seqs, 1);
        regionInd = randi([1, length(from)], addSeqsAmount, 1);
        address = rand(addSeqsAmount, 1);
        starts = from(regionInd) + floor((to(regionInd) - from(regionInd) - L) .* address);
        % load('/cs/cbio/tommy/Enhancers/Data/genome_mm9.mat');
        newSeqs = getSeqs(starts, starts + L, genome, chromoLengthUB, chrs);
        newSeqs(sum((newSeqs > 4), 2) > 0, :) = [];
        seqs = [seqs ; newSeqs];
    end
end

function seqs = getSeqs(from, to, genome, chromoLengthUB, chrs)
    chrInd = floor(from / chromoLengthUB);
    from = mod(from, chromoLengthUB);
    to = mod(to, chromoLengthUB);
    chrNames = chrs(chrInd);
    seqs = zeros(length(from), to(1) - from(1));
    for i = 1:length(from)
        fullChromosome = getfield(genome, chrNames{i});
        seqs(i, :) = nt2int(upper(fullChromosome(from(i):to(i)-1)));
    end
end


function [from, to] = merge(from, to)
    [~, ord] = sort(from);
    % sort
    from = from(ord);
    to = to(ord);
    
    M = 0;
    while M ~= size(from, 1)
        M = size(from, 1);

        % mark to merge
        toMerge = false(M, 1);
        toMerge(2:end) = to(1:end-1) > from(2:end);
        toMerge(2:end) = diff(toMerge) > 0;

        % merge
        to([toMerge(2:end);false]) = max([to([toMerge(2:end); false]), to(toMerge)], [], 2);

        % remove merged
        from(toMerge) = [];
        to(toMerge) = [];
    end
end

function [from, to] = readGenes(genes, chromoLengthUB, chrs)
    % load('/cs/cbio/tommy/Enhancers/Data/GSE29184_Bing_Ren/knownGeneXref.mm9.mat');
    fromText = {genes.ex_start}.';
    toText = {genes.ex_end}.';
    [~,chrInd] = ismember({genes.chr}, chrs);
    % unique(chrInd)
    % a = {genes.chr};
    % unique(a(chrInd ~= 0))
    fromText(chrInd == 0) = [];
    toText(chrInd == 0) = [];
    chrInd(chrInd == 0) = [];

    chrAddition = num2cell(chromoLengthUB .* chrInd).';
    from = cellfun(@str2num, fromText, 'UniformOutput', false);
    from = cellfun(@plus, from, chrAddition, 'UniformOutput', false);
    from = [from{:}].';

    to = cellfun(@str2num, toText,'UniformOutput', false);
    to = cellfun(@plus, to, chrAddition, 'UniformOutput', false);
    to = [to{:}].';
    
end