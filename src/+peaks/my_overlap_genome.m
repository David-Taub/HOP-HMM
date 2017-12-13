% generates enhancer and non enhancer datasets with overlap table
% uses ChIP-Seq from multiple tissues from MM9 genome
% genomePath = '/cs/cbio/tommy/Enhancers/Data/genome_mm9.mat';
% load(genomePath);
% knownGenesPath = '/cs/cbio/tommy/Enhancers/Data/knownGene.txt';
% [~,genesChr,~,genesStart,genesEnd] = textread(knownGenesPath,'%s%s%s%d%d%*[^\n]','headerlines',0,'delimiter','\t','bufsize',1e6);
% bedGraph = '/cs/cbio/tommy/Enhancers/Data/crgMapabilityAlign50mer100.bedGraph';
% [uniqueMapped.chr, uniqueMapped.start, uniqueMapped.end, uniqueMapped.isUnique] = textread(bedGraph, '%s%d%d%d%*[^\n]','headerlines',0,'delimiter','\t','bufsize',1e6);
% addpath('/cs/stud/boogalla/projects/CompGenetics/BaumWelch/src/data/peaks')
%  my_overlap_genome(genome, genesChr,genesStart,genesEnd, uniqueMapped)
function my_overlap_genome(genome, genesChr, genesStart, genesEnd, uniqueMapped)

    % knownGenesPath = fullfile('../data', 'knownGene.txt');
    % [~,genesChr,~,genesStart,genesEnd] = textread(knownGenesPath,'%s%s%s%d%d%*[^\n]','headerlines',0,'delimiter','\t','bufsize',1e6);
    % [uniqueMapped.chr, uniqueMapped.start, uniqueMapped.end, uniqueMapped.isUnique] = textread('crgMapabilityAlign50mer100.bedGraph','%s%d%d%d%*[^\n]','headerlines',0,'delimiter','\t','bufsize',1e6);


    global N_BIT  GENE_MARGIN_BIT  P300_BIT  K27AC_BIT  UNIQUE_BIT
    N_BIT = 1;
    GENE_MARGIN_BIT = 2;
    P300_BIT = 3;
    K27AC_BIT = 4;
    UNIQUE_BIT = 5;
    minSeqLength = 500;
    basePath = '/cs/cbio/tommy/Enhancers/Data';
    chipDirs = getChipDirs(basePath);
    p300FilePaths = getFilePathsByWord(chipDirs, 'P300', basePath);
    K27FilePaths = getFilePathsByWord(chipDirs, 'K27AC', basePath);

    chrs = fieldnames(genome);
    bitMap = buildBitMap(chrs, genome);
    genome = upperLetters(genome, chrs);
    bitMap = markN(genome, chrs, bitMap, N_BIT);
    bitMap = genesMargins(genome, bitMap, genesChr, genesStart, genesEnd, GENE_MARGIN_BIT);
    bitMap = addPeaks(bitMap, p300FilePaths, P300_BIT);
    bitMap = addPeaks(bitMap, K27FilePaths, K27AC_BIT);
    bitMap = addUniqueMapped(bitMap, uniqueMapped, UNIQUE_BIT);

    Enhancers = cell(0);
    NEnhancers = cell(0);
    for i=1:length(chrs)
        chrName = chrs{i};
        EnhancerMask = ~bitget(bitMap.(chrName),N_BIT) &...
                       ~bitget(bitMap.(chrName),GENE_MARGIN_BIT) &...
                        bitget(bitMap.(chrName),P300_BIT) &...
                        bitget(bitMap.(chrName),K27AC_BIT) &...
                        bitget(bitMap.(chrName),UNIQUE_BIT);

        NEnhancerMask = ~bitget(bitMap.(chrName),N_BIT) &...
                        ~bitget(bitMap.(chrName),GENE_MARGIN_BIT) &...
                        ~bitget(bitMap.(chrName),P300_BIT) &...
                        ~bitget(bitMap.(chrName),K27AC_BIT) &...
                         bitget(bitMap.(chrName),UNIQUE_BIT);
        % Enhancers = buildSequences(EnhancerMask, minSeqLength, chrs, Enhancers);
        % NEnhancers = buildSequences(NEnhancerMask, minSeqLength, chrs, NEnhancers);
        fprintf('%s, %f %f %f %f %f %f %f \n',  chrName,...
                    mean(bitget(bitMap.(chrName), N_BIT)),...
                    mean(bitget(bitMap.(chrName), GENE_MARGIN_BIT)),...
                    mean(bitget(bitMap.(chrName), P300_BIT)),...
                    mean(bitget(bitMap.(chrName), K27AC_BIT)),...
                    mean(bitget(bitMap.(chrName), UNIQUE_BIT)),...
                    mean(NEnhancerMask),...
                    mean(EnhancerMask));
    end
    % keyboard
    Enhancers = annotate_peaks(Enhancers);10
    NEnhancers = annotate_peaks(NEnhancers);11

    % save('-v7.3','Enhancer_map_mm9.mat','bitMap','chrs','p300','K27');
    save('Enhancers.mat','Enhancers');12
    save('NEnhancers.mat','NEnhancers');13




    % mask H3K27ac peaks


    % mark unique regions (100% unique, up to two mismatches, read length of 50)




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

end

function bitMap = buildBitMap(chrs, genome)
    fprintf('Building empty bitmap\n');
    for i=1:length(chrs),
        chrName=chrs{i};
        len(i)=length(genome.(chrName));
        bitMap.(chrName)=uint8(zeros(1,len(i)));
    end
end

function genome = upperLetters(genome, chrs)
    fprintf('Upper casing chromosome names \n');
    for i=1:length(chrs),
        chrName =chrs{i};
        genome.(chrName)=upper(genome.(chrName));
    end
end

function bitMap = markN(genome, chrs, bitMap, bitNum)
    fprintf('N marking\n');
    % mask out Ns
    for i=1:length(chrs),
        chrName = chrs{i};
        bitMap.(chrName)(genome.(chrName)=='N') = bitNum;

    end
end

function bitMap = genesMargins(genome, bitMap, genesChr, genesStart, genesEnd, bitNum)
    fprintf('Gene with margins marking\n');
    % mask out genes +/- 15000
    % [~,chrName,~,genesStart,genesEnd = textread('knownGene.txt','%s%s%s%d%d%*[^\n]','headerlines',0,'delimiter','\t','bufsize',1e6);
    for i=1:length(genesChr),
        regionStart = max(1,genesStart(i)+1-15000);
        regionEnd = min(genesEnd(i)+15000,length(genome.(genesChr{i})));
        bitMap.(genesChr{i})(regionStart:regionEnd) = bitset(bitMap.(genesChr{i})(regionStart:regionEnd), bitNum);
    end
end

function chipDirs = getChipDirs(basePath)
    fprintf('get all ChIP dirs\n');
    chipDirs = dir(fullfile(basePath, 'G*'));
    for i=length(chipDirs):-1:1
        if ~chipDirs(i).isdir
            chipDirs(i)=[];
        end;
    end
end

function filePaths = getFilePathsByWord(chipDirs, word, basePath)
    fprintf('get all ChIP files %s\n', word);
    filePaths={};
    for i=1:length(chipDirs)
        dirName=chipDirs(i).name;
        matFiles = dir(fullfile(basePath, dirName, '*.mat'));
        for j=1:length(matFiles),
            if isempty(strfind(upper(matFiles(j).name),word))
                continue;
            end;
            filePaths{end+1}= fullfile(basePath, dirName, matFiles(j).name);
        end
    end
    fprintf('found %d files\n', length(filePaths));
end

function bitMap = addPeaks(bitMap, filePaths, bitNum)

    fprintf('peaks marking, %d files \n', length(filePaths));
    for i = 1:length(filePaths)
        A = load(filePaths{i});
        fprintf('%s %d\n', filePaths{i}, length(A.S));
        for j = 1:length(A.S)
            peak = A.S{j};
            bitMap.(peak.chr)(peak.from:peak.to) = bitset(bitMap.(peak.chr)(peak.from:peak.to),bitNum);
        end
    end
end
function bitMap = addUniqueMapped(bitMap, uniqueMapped, bitNum)
    fprintf('mark unique loci\n');
    % [uniqueMapped.chr, uniqueMapped.start, uniqueMapped.end, uniqueMapped.isUnique] = textread('crgMapabilityAlign50mer100.bedGraph','%s%d%d%d%*[^\n]','headerlines',0,'delimiter','\t','bufsize',1e6);
    indices=find(uniqueMapped.isUnique==1);
    for i=1:length(indices),
        j=indices(i);
        regionStart = uniqueMapped.start(j)+1;
        regionEnd = uniqueMapped.end(j); % uniquely-mapped block
        bitMap.(uniqueMapped.chr{j})(regionStart:regionEnd) = bitset(bitMap.(uniqueMapped.chr{j})(regionStart:regionEnd),bitNum);
    end
end

function sequences = buildSequences(mask, genome, minSeqLength, chrs, sequences)
    fprintf('building sequences structs\n');
    sequences = {};
    for i=1:length(chrs)
        chrName=chrs{i};
        mask = 2*mask-1;

        startsTmp = conv(mask,[1 1 -1]);
        startsTmp = startsTmp(2:end-1);
        starts = find(startsTmp==3);
        if startsTmp(1)==2
            starts = [1 starts];
        end

        endsTmp = conv(mask,[1 -1 -1]);
        endsTmp = endsTmp(2:end-1);
        ends = find(endsTmp == -3');

        for j=1:length(starts)
            if ends(j)-starts(j)<minSeqLength
                continue;
            end
            seq.chr = chrName;
            seq.from = starts(j);
            seq.to = ends(j);
            seq.seq = upper(genome.(chrName)(starts(j):ends(j)))
            sequences{end+1} = seq;
        end
    end
end
