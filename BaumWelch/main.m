function main()
    close all;
    enTrain = readSeq('Data/Enhancers.train.seq');
    nenTrain = readSeq('Data/NEnhancers.train.seq');
    enTest = readSeq('Data/Enhancers.test.seq');
    nenTest = readSeq('Data/NEnhancers.test.seq');

    s = 0.0005;
    for order = 1 : 4
        Een = getEFromSeqs(enTrain, order);
        Enen = getEFromSeqs(nenTrain, order);
        logLikes1 = getLogLikes(Een, enTrain, order);
        logLikes2 = getLogLikes(Enen, enTrain, order);
        logLikes3 = getLogLikes(Een, nenTrain, order);
        logLikes4 = getLogLikes(Enen, nenTrain, order);

        ratioEn = logLikes1 ./ logLikes2;
        ratioNen = logLikes3 ./ logLikes4;


        likeMin = min([ratioEn; ratioNen]);
        likeMax = max([ratioEn; ratioNen]);
        Ls = [likeMin : s : likeMax];
        subplot(2,2,order); 
        plot(Ls, histc(ratioEn, Ls));
        hold on
        plot(Ls, histc(ratioNen, Ls));
        legend('Ratio for en', 'Ratio for Non en');
        title(sprintf('Order %d', order))
        hold off
    end
end

function E = getEFromSeqs(seqs, order)
    N = length(seqs);
    matSize = [4 * ones(1, order), 1];
    E = zeros(matSize);
    Es = zeros(prod(4 ^ order), N);
    for i=1:N
        fprintf('Getting emission matrix %d / %d\r', i, N)
        Ecur = getEFromSeq(seqs{i}, order);
        Es(:, i) = Ecur(:);
    end
    fprintf('\nDone.\n');
    E(:) = mean(Es, 2);
end

function logLikes = getLogLikes(E, seqs, order)
    N = length(seqs);
    logLikes = zeros(N,1);
    for i=1:N
        fprintf('Getting log likelihood %d / %d\r', i, N)
        logLikes(i) = getLogLikeFromSeq(seqs{i}, E, order);
    end
    fprintf('\nDone.\n');
end

% reads .seq format file (Tommy's format) and returns the sequence as numbers
function seq = readSeq(filePath)
    fprintf('Reading %s...', filePath);
    [acc, seq] = textread(filePath, '%s%s%*[^\n]','delimiter','\t','bufsize',20000);
    for i=1:length(seq)
        seq{i}=upper(seq{i}); 
    end;
    N = length(seq);

    len = zeros(1,N); for i=1:N, len(i)=length(seq{i}); end;
    for i=1:length(seq), seq{i}=seq{i}(ceil(len(i)/2) + [-249:250]); end;

    % map sequences to dinucleotide indices
    for i=1:N,
        seq{i} = nt2int(seq{i});
    end
    fprintf(' done.\n');
end

% very similar to sub2ind, but receives the subscripts as matrix
% matSize - 1 x k
% subscripts - k x n
% indices - 1 x n
function indices = mySub2ind(matSize, subscripts)
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
    indices = mySub2ind(matSize, k);
end

% seq - 1 x n
% E - 4 x ... x 4 (order times)
function E = getEFromSeq(seq, order)
    ambient = 10 ^ -6;
    matSize = 4 * ones(1, order);
    % indices - 1 x n
    indices = getIndeices1D(seq, order);
    h = histc(indices, 1 : 4 ^ order);
    h = h + ambient;
    E = reshape(h, [matSize, 1]);
    E = bsxfun(@times, E, 1 ./ sum(E, order));
end

% seq - 1 x n
% E - 4 x ... x 4 (order times)
% logLike - number
function logLike = getLogLikeFromSeq(seq, E, order)
    indices = getIndeices1D(seq, order);
    logLike = sum(sum(log(E(indices)), 1), 2);
    % logLike
    logLike = exp(logLike / length(seq));
end