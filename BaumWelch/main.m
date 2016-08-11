function main()
    close all;
    L = 700;
    posSeqsTrain = readSeq('Enhancers.train.seq', L);
    negSeqsTrain = readSeq('NEnhancers.train.seq', L);
    posSeqsTest = readSeq('Enhancers.test.seq', L);
    negSeqsTest = readSeq('NEnhancers.test.seq', L);
    bestTrainError = inf;
    bestOrder = 0;
    bestAlpha = 0;
    bestEen = 0;
    bestEnen = 0;
    for order = 1:8
        
        [posE, negE] = trainMarkov(posSeqsTrain, negSeqsTrain, order);
        
        thresholds = 0.98 : 0.00005 : 1.01;
        % thresholds = 1.0004;
        trainErrs = classify(posE, negE, posSeqsTrain, negSeqsTrain, thresholds);
        [trainErr, i] = min(trainErrs);
        threshold = thresholds(i);
        testErr = classify(posE, negE, posSeqsTest, negSeqsTest, threshold);
        [order, threshold, trainErr, testErr]
    end
end

function [posE, negE] = trainMarkov(posSeqs, negSeqs, order)
    posE = getEFromSeqs(posSeqs, order);
    negE = getEFromSeqs(negSeqs, order);
end

function err = classify(posE, negE, posSeqs, negSeqs, thresholds)
    likePosIsPos = getLogLikes(posE, posSeqs); %high
    likePosIsNeg = getLogLikes(negE, posSeqs); %low
    likeNegIsPos = getLogLikes(posE, negSeqs); %low
    likeNegIsNeg = getLogLikes(negE, negSeqs); %high

    ratioPos = likePosIsPos ./ likePosIsNeg; %high / low = high
    ratioNeg = likeNegIsPos ./ likeNegIsNeg; %low
    % m = min([likePosIsPos ; likePosIsNeg ; likeNegIsPos ; likeNegIsNeg])
    % M = max([likePosIsPos ; likePosIsNeg ; likeNegIsPos ; likeNegIsNeg])
    % Ls = m : 1: M;
    % hold on
    % plot(Ls, histc(likePosIsPos, Ls));
    % plot(Ls, histc(likePosIsNeg, Ls));
    % plot(Ls, histc(likeNegIsPos, Ls));
    % plot(Ls, histc(likeNegIsNeg, Ls));
    % legend('pp', 'pn', 'np', 'nn')
    % hold off;
    
    err = zeros(size(thresholds));
    for i = 1:length(thresholds);
        err(i) = sum(ratioPos > thresholds(i)) + ...
                 sum(ratioNeg < thresholds(i));
    end
    err = err ./ (length(ratioPos) + length(ratioNeg));
end

function E = getEFromSeqs(seqs, order)
    ambient = 10 ^ -6;
    N = length(seqs);
    matSize = [4 * ones(1, order), 1];
    E = zeros(matSize);
    Es = zeros(prod(4 ^ order), N);
    for i = 1 : N
        % fprintf('Getting emission matrix %d / %d\r', i, N)
        indices = getIndeices1D(seqs{i}, order);
        h = histc(indices, 1 : 4 ^ order);
        Ecur = reshape(h, [matSize, 1]);
        E = E + Ecur;
    end
    h = h + ambient;
    E = bsxfun(@times, E, 1 ./ sum(E, order));
    % fprintf('\nDone.\n');
end

function logLikes = getLogLikes(E, seqs)
    N = length(seqs);
    logLikes = zeros(N,1);
    for i=1:N
        % fprintf('Getting log likelihood %d / %d\r', i, N)
        logLikes(i) = getLogLikeFromSeq(seqs{i}, E);
    end
    % fprintf('\nDone.\n');
end

% reads .seq format file (Tommy's format) and returns the sequence as numbers
function seq = readSeq(filePath, L)
    % fprintf('Reading %s...', filePath);
    [acc, seq] = textread(filePath, '%s%s%*[^\n]','delimiter','\t','bufsize',20000);
    for i=1:length(seq)
        seq{i}=upper(seq{i}); 
    end;
    N = length(seq);
    len = zeros(1,N); for i=1:N, len(i)=length(seq{i}); end;
    seq(len < L) = [];
    len(len < L) = [];
    N = length(seq);
    for i=1:length(seq), seq{i}=seq{i}(ceil(len(i)/2) + [-floor(L/2) + 1: floor(L/2)]); end;

    % map sequences to dinucleotide indices
    for i=1:N,
        seq{i} = nt2int(seq{i});
    end
    % fprintf(' done.\n');
end

% very similar to sub2ind, but receives the subscripts as matrix
% matSize - 1 x k
% subscripts - k x n
% indices - 1 x n
function indices = matSub2ind(matSize, subscripts)
    subtractedSub = subscripts.';
    subtractedSub(2:end, :) = subtractedSub(2:end, :) - 1;
    % subtractedSub
    cumMatSize = cumprod([1, matSize]);
    cumMatSize = cumMatSize(1, 1:end - 1);
    indices = cumMatSize * subtractedSub;
end


% seq - 1 x n
% indices - 1 x n (numbers from 1 to order)
function indices = getIndeices1D(seq, order)
    N = length(seq);
    k = zeros(N - order + 1, order);
    for i = 1 : order
        seq(i : end - order + i);
        k(:,i) = seq(i : end - order + i);
    end
    matSize = 4 * ones(1, order);
    indices = matSub2ind(matSize, k);
end

% seq - 1 x n
% E - 4 x ... x 4 (order times)
% logLike - number
function logLike = getLogLikeFromSeq(seq, E)
    order = matDim(E);
    indices = getIndeices1D(seq, order);
    logLike = sum(sum(log(E(indices)), 1), 2);
    % logLike
    logLike = logLike;
end
