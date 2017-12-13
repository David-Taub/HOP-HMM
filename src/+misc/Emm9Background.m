% misc.Emm9Background();

% E - 4 x 4 x 4
function E = Emm9Background()
    n = 4;
    L = 500;
    order = 3;
    [~, x] = loadTommySeqs(L);
    x = x(:);
    x = x(1:end-mod(length(x), order));
    x = reshape(x, [length(x)./order, order]);
    x = x * [1; 10; 100];
    xx = unique(x);
    xx = sort(xx);          % sorted input aligns with temp (lowest to highest)
    t = zeros(size(xx));
    for i = 1:length(xx)
        t(i) = sum(x == xx(i));
    end
    ESize = ones(1,order) * n;
    preE = reshape(t, ESize);
    E = preE ./ repmat(sum(preE, order), [ones(1,order-1), n]);
    E = permute(E, [order+1, 1:order]);
    save('data/temp/mm9NonEnhE.mat','E');
end
% reads .seq format file (Tommy's format) and returns the sequence as numbers
function out = readSeq(filePath, L)
    [~, seqsCells] = textread(filePath, '%s%s%*[^\n]','delimiter','\t','bufsize',20000);
    out = regularSeqs(seqsCells, L);
    % save(matPath, 'out');
end

% L - sequence lengths
% n - number of classes to get for the positive sequences, where class
% means unique overlap between tissues, and if n is 1:3 then we take
% the sequences of the three most frequent class
function [posSeqs, negSeqs] = loadTommySeqs(L)
    negSeqsTrain = readSeq('/cs/cbio/tommy/Enhancers/Data/NEnhancers.train.seq', L);
    posSeqsTrain = readSeq('/cs/cbio/tommy/Enhancers/Data/Enhancers.train.seq', L);
    negSeqsTest = readSeq('/cs/cbio/tommy/Enhancers/Data/NEnhancers.test.seq', L);
    posSeqsTest = readSeq('/cs/cbio/tommy/Enhancers/Data/Enhancers.test.seq', L);

    posSeqs = [posSeqsTest; posSeqsTrain];
    negSeqs = [negSeqsTest; negSeqsTrain];

    N = min(size(negSeqs, 1), size(posSeqs, 1));
    posSeqs = posSeqs(1:N, :);
    negSeqs = negSeqs(1:N, :);
    % shuffle
    negSeqs = negSeqs(randperm(size(negSeqs, 1)), :);
    posSeqs = posSeqs(randperm(size(posSeqs, 1)), :);
end


% remove short seqs
% ACGT -> 1234
function out = regularSeqs(seqsCells, L)
    out = zeros(length(seqsCells), L);
    for i=1:length(seqsCells)
        seqsCells{i}=upper(seqsCells{i});
        if length(seqsCells{i}) < L
            continue;
        end
        seqsCells{i}=seqsCells{i}(ceil(length(seqsCells{i})/2) + [-floor(L/2) + 1: floor(L/2)]);
        out(i, :) = nt2int(seqsCells{i});
    end
    out( ~any(out,2), : ) = [];  %remove zero rows
    % in some of the data we have NNN which can be any nucleotide
    out(out==15) = 1;
end