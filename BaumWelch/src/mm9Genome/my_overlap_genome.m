function [] = my_overlap_genome(genome)

    % load('genome_mm9.mat');
    chrs = fieldnames(genome);

    fprintf('Upper casing chromosome names \n');
    for i=1:length(chrs),
        c=chrs{i};
        genome.(c)=upper(genome.(c));
    end

    fprintf('Building bitmap \n');
    for i=1:length(chrs),
        c=chrs{i}; 
        len(i)=length(genome.(c)); 
        bitMap.(c)=uint8(zeros(1,len(i)));
    end
    clear i c;


    fprintf('N marking\n');
    % mask out Ns
    for i=1:length(chrs),
        c = chrs{i}; 
        bitMap.(c)(genome.(c)=='N') = 2^0;
        fprintf('%s, %f %f %f %f %f\n',  c,...
                mean(bitget(bitMap.(c), 1)),...
                mean(bitget(bitMap.(c), 2)),...
                mean(bitget(bitMap.(c), 3)),...
                mean(bitget(bitMap.(c), 4)),...
                mean(bitget(bitMap.(c), 5)));
    end
    clear i c;

    fprintf('Gene with margins marking\n');
    % mask out genes +/- 15000?
    [~,chrName,~,genesStart,genesEnd] = textread('knownGene.txt','%s%s%s%d%d%*[^\n]','headerlines',0,'delimiter','\t','bufsize',1e6);
    for i=1:length(chrName),
        regionStart = max(1,genesStart(i)+1-15000); 
        regionEnd = min(genesEnd(i)+15000,length(genome.(chrName{i})));
        bitMap.(chrName{i})(regionStart:regionEnd) = bitset(bitMap.(chrName{i})(regionStart:regionEnd),2);
    end
    clear chrName genesStart genesEnd i;


    % get all ChIP files
    dirs = dir('G*'); 
    ld=length(dirs); 
    for i=ld:-1:1
        if ~isdir(dirs(i).name)
            dirs(i)=[];
        end;
        ld=length(dirs); 
    end;


    % mask p300 peaks
    p300={};
    for i=1:ld
        d=dirs(i).name; 
        f=dir([d '/*.mat']);
        for j=1:length(f),
            if isempty(strfind(upper(f(j).name),'P300'))
                continue;
            end;
            p300{end+1}=[d '/' f(j).name];
        end;
    end;

    fprintf('p300 peaks marking, %d files \n', length(p300));
    for i = 1:length(p300)
        A = load(p300{i});
        fprintf('%s %d\n', p300{i}, length(A.S));
        for j = 1:length(A.S)
            peak = A.S{i};
            bitMap.(peak.chr)(peak.from:peak.to) = bitset(bitMap.(peak.chr)(peak.from:peak.to),3);
        end
    end
    clear i j d f A peak


    % mask H3K27ac peaks
    K27={};
    for i=1:ld,
        d=dirs(i).name; 
        f=dir([d '/*.mat']);
        for j=1:length(f),
            if isempty(strfind(upper(f(j).name),'K27AC'))
                continue;
            end;
            K27{end+1}=[d '/' f(j).name];
        end;
    end;

    fprintf('H3K27ac peaks marking, %d files \n', length(K27));
    for i=1:length(K27),
        A=load(K27{i});
        if isfield(A, 'S')
            fprintf('%s %d \n', K27{i}, length(A.S));
            for i=1:length(A.S),
                peak = A.S{i}; 
                bitMap.(peak.chr)(peak.from:peak.to) = bitset(bitMap.(peak.chr)(peak.from:peak.to),4);
            end
        end
    end
    clear i j d f A peak

    % mark unique regions (100% unique, up to two mismatches, read length of 50)
    [c,f,t,m] = textread('crgMapabilityAlign50mer100.bedGraph','%s%d%d%d%*[^\n]','headerlines',0,'delimiter','\t','bufsize',1e6);
    I=find(m==1);
    for i=1:length(I),
        j=I(i); 
        regionStart = f(j)+1; 
        regionEnd = t(j); % uniquely-mapped block
        bitMap.(c{j})(regionStart:regionEnd) = bitset(bitMap.(c{j})(regionStart:regionEnd),5);
    end
    clear c f t m i j I;



    % mask H3K4me1 peaks


    % bitmap
    % 1 = N
    % 2 = within gene (up to 15K, including introns)
    % 3 = p300 peak (within any condition)
    % 4 = H3K27ac peak (within any condition)
    % 5 = Unique region (100% unique, >=2 mismatches, read length of 50)
    
    %                   N           gene15K     p300        H3K27ac     Unique
    %
    % chr1,             0.028997    0.497680    0.000021    0.191179    0.812995
    % chr2,             0.018465    0.606896    0.000000    0.235026    0.833892
    % chr3,             0.020087    0.450341    0.000000    0.178025    0.808633
    % chr4,             0.024052    0.549327    0.000000    0.212842    0.796101
    % chr5,             0.031573    0.573740    0.000060    0.213289    0.803958
    % chr6,             0.021402    0.593912    0.000000    0.207700    0.810489
    % chr7,             0.069801    0.591607    0.000014    0.196058    0.706706
    % chr8,             0.052696    0.499503    0.000000    0.210511    0.806093
    % chr9,             0.027047    0.598585    0.000003    0.243316    0.838326
    % chr10,            0.024197    0.537986    0.000000    0.216998    0.825983
    % chr11,            0.025445    0.636968    0.000000    0.287201    0.857547
    % chr12,            0.031324    0.479885    0.000000    0.196779    0.772933
    % chr13,            0.032535    0.508708    0.000000    0.216152    0.789561
    % chr14,            0.028432    0.518395    0.000000    0.168887    0.762074
    % chr15,            0.029518    0.518560    0.000000    0.216736    0.828513
    % chr16,            0.033709    0.517970    0.000022    0.196625    0.824768
    % chr17,            0.035419    0.571779    0.000000    0.224281    0.777536
    % chr18,            0.034944    0.484831    0.000000    0.211895    0.826140
    % chr19,            0.052169    0.631623    0.000000    0.252758    0.814134
    % chrM,             0.000000    1.000000    0.000000    0.000000    0.480888
    % chrX,             0.027419    0.409380    0.000000    0.059770    0.652134
    % chrY,             0.830055    0.083318    0.000000    0.001791    0.068559

    % chr1_random,      0.180966    0.479207    0.000000    0.048998    0.198408
    % chr3_random,      0.000000    0.000000    0.000000    0.000000    0.001456
    % chr4_random,      0.329850    0.313418    0.000000    0.035045    0.357105
    % chr5_random,      0.139919    0.295265    0.000000    0.000000    0.000327
    % chr7_random,      0.137935    0.289065    0.000000    0.000000    0.006596
    % chr8_random,      0.176555    0.529157    0.000000    0.001647    0.055769
    % chr9_random,      0.148508    0.533007    0.000000    0.063362    0.617034
    % chr13_random,     0.160990    0.826422    0.000000    0.008546    0.095064
    % chr16_random,     0.000000    0.000000    0.000000    0.000000    0.000000
    % chr17_random,     0.159049    0.497777    0.000000    0.000000    0.000463
    % chrX_random,      0.292957    0.353627    0.000000    0.004225    0.125434
    % chrY_random,      0.086056    0.072674    0.000000    0.000000    0.012076
    % chrUn_random,     0.458238    0.188948    0.000000    0.000832    0.102224
% GSE13845_Forebrain_Midbrain_Limb/GSM348064_forebrain_p300_peaks.mat 2453
% GSE13845_Forebrain_Midbrain_Limb/GSM348065_midbrain_p300_peaks.mat 561
% GSE13845_Forebrain_Midbrain_Limb/GSM348066_limb_p300_peaks.mat 2105
% GSE32587_Heart/GSM807737_mouse_p2_heart_accbp-p300_peak.mat 6565
% GSE36027_GSE31039_ENCODE/GSM918747_mm9_wgEncodeLicrTfbsHeartP300MAdult8wksC57bl6StdPk.mat 39447
% GSE36027_GSE31039_ENCODE/GSM918750_mm9_wgEncodeLicrTfbsEsb4P300ME0C57bl6StdPk.mat 23654
% GSE37151_catalog/GSM1135064_mouse_e11.5_face_ep300-flag.peaks.mat 3589
% GSE37151_catalog/GSM1135066_mouse_ESC_ep300-flag.peaks.mat 10356
% GSE42881_Brain/GSM1052708_mouse_e11.5_forebrain_p300.peaks.mat 4425
% GSE42881_Brain/GSM1052709_mouse_p0_cortex_p300.peaks.mat 1132
% GSE49413_Craniofacial/GSM1199037_mm9_p300_craniofacial_e11.5_peaks.mat 7282


    % for example - this could serve as a negative control for enhancer regions
    % no N's, not near/in genes, no known enhancer marks, and unique (34% of genome?)
    for i=1:length(chrs),
            c = chrs{i}; 
        fprintf('%s, %f %f %f %f %f\n',  c,...
                mean(bitget(bitMap.(c), 1)),...
                mean(bitget(bitMap.(c), 2)),...
                mean(bitget(bitMap.(c), 3)),...
                mean(bitget(bitMap.(c), 4)),...
                mean(bitget(bitMap.(c), 5)));
    end

    mean(   ~bitget(bitMap.chr1, 1) &...
            ~bitget(bitMap.chr1, 2) &...
            ~bitget(bitMap.chr1, 3) &...
            ~bitget(bitMap.chr1, 4) &...
             bitget(bitMap.chr1, 5))

    % Enhancers
    % no N's, not near/in genes, known enhancer marks, unique, and at least 500bp long
    segs = []; Enhancers = cell(0);
    for i=1:length(chrs),
        c=chrs{i};
        A = ~bitget(bitMap.(c),1) &...
            ~bitget(bitMap.(c),2) &...
             bitget(bitMap.(c),3) &...
             bitget(bitMap.(c),4) &...
             bitget(bitMap.(c),5);
        regionStart=2*A-1;
        AA = conv(regionStart,[1 1 -1]); AA = AA(2:end-1); st = find(AA==3); if AA(1)==2, st = [1 st]; end;
        AA = conv(regionStart,[1 -1 -1]); AA = AA(2:end-1); en = find(AA==-3'); 
        segs.(c) = [st;en]';

        for j=1:length(st),
    	if en(j)-st(j)<500, continue; end;
    	clear s; s.chr=c; s.from=st(j); s.to=en(j);
    	Enhancers{end+1}=s;
        end
    end

    % Non-Enhancers
    % no N's, not near/in genes, no known enhancer marks, unique, and at least 500bp long
    segs = []; NEnhancers = cell(0);
    for i=1:length(chrs),
        c=chrs{i};
        A = ~bitget(bitMap.(c),1) & ~bitget(bitMap.(c),2) & ~bitget(bitMap.(c),3) & ~bitget(bitMap.(c),4) & bitget(bitMap.(c),5);
        regionStart=2*A-1;
        AA = conv(regionStart,[1 1 -1]); AA = AA(2:end-1); st = find(AA==3); if AA(1)==2, st = [1 st]; end;
        AA = conv(regionStart,[1 -1 -1]); AA = AA(2:end-1); en = find(AA==-3'); 
        segs.(c) = [st;en]';

        for j=1:length(st),
    	if en(j)-st(j)<500, continue; end;
    	clear s; s.chr=c; s.from=st(j); s.to=en(j);
    	NEnhancers{end+1}=s;
        end
    end
    clear segs i c A regionStart j st en

    Enhancers = annotate_peaks(Enhancers);
    NEnhancers = annotate_peaks(NEnhancers);

    % load('genome_mm9.mat');
    for i=1:length(Enhancers), s=Enhancers{i}; c=s.chr; f=s.from; t=s.to; Enhancers{i}.seq=upper(genome.(c)(f:t)); end
    for i=1:length(NEnhancers), s=NEnhancers{i}; c=s.chr; f=s.from; t=s.to; NEnhancers{i}.seq=upper(genome.(c)(f:t)); end
    clear i s c f t seq

    save('-v7.3','Enhancer_map_mm9.mat','bitMap','chrs','p300','K27');
    save('Enhancer_map_mm9_Enhancers.mat','Enhancers');
    save('Enhancer_map_mm9_NEnhancers.mat','NEnhancers');
end