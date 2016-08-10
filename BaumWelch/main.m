function main()
    enTrain = readSeq('Data/Enhancers.train.seq');
    nenTrain = readSeq('Data/NEnhancers.train.seq');
    enTest = readSeq('Data/Enhancers.test.seq');
    nenTest = readSeq('Data/NEnhancers.test.seq');
    getE(2, enTrain);   
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
function ind = sub2ind(seq)
    for i=1:n
% matSize - 1 x k
% subscripts - k x n
% indices - 1 x n
function indices = mySub2ind(matSize, subscripts)
    subtractedSub = subscripts;
    subtractedSub(2:end, :) = subtractedSub(2:end, :) - 1;
    cumMatSize = cumprod([1, matSize]);
    cumMatSize = cumMatSize(1, end -1);
    indices = cumMatSize * subtractedSub;
end


function E = getE(order, seq)
    N = length(seq);
    k = zeros(N - order + 1, order);
    for i = 1 : order
        k(:,i) = seq(i : end - order + i);
    end

    % map sequences to dinucleotide indices
    matSize = 4 * ones(1, order);
    indices = mySub2ind(matSize, k);
    E = reshape(histc(indices, 1 : 4 ^ order), matSize);


        % pentnuc{i} = sub2ind([4,4,4,4,4,4,4],t(7:end),t(6:end-1),t(5:end-2),t(4:end-3),t(3:end-4),t(2:end-5),t(1:end-6));

    for i=1:4,   mono{i} = int2nt(dec2base(i-1,4,1) - 'a' + 50); end;
    for i=1:4^2, di{i}   = int2nt(dec2base(i-1,4,2) - 'a' + 50); end;
    for i=1:4^3, tri{i}  = int2nt(dec2base(i-1,4,3) - 'a' + 50); end;
    for i=1:4^4, quad{i} = int2nt(dec2base(i-1,4,4) - 'a' + 50); end;
    for i=1:4^5, pent{i} = int2nt(dec2base(i-1,4,5) - 'a' + 50); end;
    for i=1:4^6, hexa{i} = int2nt(dec2base(i-1,4,6) - 'a' + 50); end;
    for i=1:4^7, pent{i} = int2nt(dec2base(i-1,4,7) - 'a' + 50); end;

    alpha = 0.01;
    MONOPCT = zeros(N,4);   for i=1:N, MONOPCT(i,:) = (alpha + histc(mononuc{i},1:4))  / ((4   * alpha) + len(i)); end;
    DIPCT   = zeros(N,4^2); for i=1:N, DIPCT(i,:)   = (alpha + histc(dinuc{i},1:4^2))  / ((4^2 * alpha) + (len(i)-1)); end;
    TRIPCT  = zeros(N,4^3); for i=1:N, TRIPCT(i,:)  = (alpha + histc(trinuc{i},1:4^3)) / ((4^3 * alpha) + (len(i)-2)); end;
    QUADPCT = zeros(N,4^4); for i=1:N, QUADPCT(i,:) = (alpha + histc(quadnuc{i},1:4^4)) / ((4^4 * alpha) + (len(i)-3)); end;
    PENTAPCT= zeros(N,4^5); for i=1:N, PENTAPCT(i,:)= (alpha + histc(pentanuc{i},1:4^5)) / ((4^5 * alpha) + (len(i)-4)); end;
    HEXAPCT = zeros(N,4^6); for i=1:N, HEXAPCT(i,:) = (alpha + histc(hexanuc{i},1:4^6)) / ((4^6 * alpha) + (len(i)-5)); end;
    % PENTPCT = zeros(N,4^7); for i=1:N, PENTPCT(i,:) = (alpha + histc(pentnuc{i},1:4^7)) / ((4^7 * alpha) + (len(i)-6)); end;

X1=MONOPCT; X2=DIPCT; X3=TRIPCT; X4=QUADPCT; X5=PENTAPCT; X6=HEXAPCT; Y=pos;


end