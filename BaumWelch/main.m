function main()
    enTrain = readSeq('Data/Enhancers.train.seq');
    nenTrain = readSeq('Data/NEnhancers.train.seq');
    enTest = readSeq('Data/Enhancers.test.seq');
    nenTest = readSeq('Data/NEnhancers.test.seq');
    getE(enTrain{1}, 3);
end

% reads .seq format file (Tommy's format) and returns the sequence as numbers
function seq = readSeq(filePath)
    [acc, seq] = textread(filePath, '%s%s%*[^\n]','delimiter','\t','bufsize',20000);
    for i=1:length(seq)
        seq{i}=upper(seq{i}); 
    end;
    N = length(seq);
    len = zeros(1,N); for i=1:N, len(i)=length(seq{i}); end;
    % try nt2int(seq{i})
    for i=1:length(seq), seq{i}=seq{i}(ceil(len(i)/2) + [-249:250]); end;
    len = zeros(1,N); for i=1:N, len(i)=length(seq{i}); end;

    % map sequences to dinucleotide indices
    for i=1:N,
        seq{i} = nt2int(seq{i});
    end
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


function E = getE(seq, order)
    N = length(seq);
    k = zeros(N - order + 1, order);
    for i = 1 : order
        seq(i : end - order + i);
        k(:,i) = seq(i : end - order + i);
    end
    matSize = 4 * ones(1, order);
    indices = mySub2ind(matSize, k);
    h = histc(indices, 1 : 4 ^ order);
    E = reshape(h, matSize);
    E = bsxfun(@times, E, 1 ./ sum(E, order));
end