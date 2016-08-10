function main()
    enTrain = readSeq('Data/Enhancers.train.seq');
    nenTrain = readSeq('Data/NEnhancers.train.seq');
    enTest = readSeq('Data/Enhancers.test.seq');
    nenTest = readSeq('Data/NEnhancers.test.seq');
    
end
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
