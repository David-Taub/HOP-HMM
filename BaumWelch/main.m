function main()
end
function getSeq()
    [accPtrain, seqPtrain] = textread('Data/Enhancers.train.seq', '%s%s%*[^\n]','delimiter','\t','bufsize',20000);
    [accNtrain, seqNtrain] = textread('Data/NEnhancers.train.seq','%s%s%*[^\n]','delimiter','\t','bufsize',20000);
    [accPtest, seqPtest]   = textread('Data/Enhancers.test.seq',  '%s%s%*[^\n]','delimiter','\t','bufsize',20000);
    [accNtest, seqNtest]   = textread('Data/NEnhancers.test.seq', '%s%s%*[^\n]','delimiter','\t','bufsize',20000);

    accP = [accPtrain; accPtest]; 
    seqP = [seqPtrain; seqPtest];
    accN = [accNtrain; accNtest]; 
    seqN = [seqNtrain; seqNtest];

    pos = [true(length(accP), 1); false(length(accN), 1)];
    neg = ~pos;
    acc = [accP; accN];
    seq = [seqP; seqN];
    train = [true(length(accPtrain),1); false(length(accPtest),1); true(length(accNtrain),1);false(length(accNtest),1)]; test=~train;

    clear accP seqP accN seqN I accPtrain seqPtrain accNtrain seqNtrain accPtest seqPtest accNtest seqNtest
    N = length(seq);
    for i=1:length(seq)
        seq{i}=upper(seq{i}); 
    end;
    len = zeros(1,N); for i=1:N, len(i)=length(seq{i}); end;
    % try nt2int(seq{i})
    for i=1:length(seq), seq{i}=seq{i}(ceil(len(i)/2) + [-249:250]); end;
    len = zeros(1,N); for i=1:N, len(i)=length(seq{i}); end;

    % map sequences to dinucleotide indices
    for i=1:N,
        t = nt2int(seq{i});
        mononuc{i} = t;
        dinuc{i} = sub2ind([4,4],t(2:end),t(1:end-1));
        trinuc{i} = sub2ind([4,4,4],t(3:end),t(2:end-1),t(1:end-2));
        quadnuc{i} = sub2ind([4,4,4,4],t(4:end),t(3:end-1),t(2:end-2),t(1:end-3));
        pentanuc{i} = sub2ind([4,4,4,4,4],t(5:end),t(4:end-1),t(3:end-2),t(2:end-3),t(1:end-4));
        hexanuc{i} = sub2ind([4,4,4,4,4,4],t(6:end),t(5:end-1),t(4:end-2),t(3:end-3),t(2:end-4),t(1:end-5));
        % pentnuc{i} = sub2ind([4,4,4,4,4,4,4],t(7:end),t(6:end-1),t(5:end-2),t(4:end-3),t(3:end-4),t(2:end-5),t(1:end-6));
    end

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

    svm2 = fitcsvm(X2(train,:),Y(train),'KernelFunction','linear','Standardize',true,'KernelScale','auto');
    [label2,score2] = predict(svm2,X2(test,:)); 100*mean(label2==Y(test)),;
    [tx2,ty2,~,auc2] = perfcurve(pos(test)',score2(:,2),true); plot(tx2,ty2,'LineWidth',2);title(100*auc2);
end
end