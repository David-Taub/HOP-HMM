function [] = overlap_genome()

% load('genome_mm9.mat');
chrs = fieldnames(genome);

for i=1:length(chrs),
    c=chrs{i};
    genome.(c)=upper(genome.(c));
end

for i=1:length(chrs),
    c=chrs{i}; 
    len(i)=length(genome.(c)); 
    bitMap.(c)=uint8(zeros(1,len(i)));
end
clear i c;

% mask out Ns
for i=1:length(chrs),
    c = chrs{i}; 
    bitMap.(c)(genome.(c)=='N') = 2^0;
end
clear i c;

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

for i = 1:length(p300)
    A = load(p300{i});
    for j = 1:length(A.S)
        peak = A.S{i};
        bitMap.(peak.chr)(peak.from:peak.to) = bitset(bitMap.(peak.chr)(peak.from:peak.to),3);
    end
end
clear i j d f A peak


% mask H3K27ac peaks
K27={};
for i=1:ld,
    d=dirs(i).name; f=dir([d '/*.mat']);
    for j=1:length(f),
        if isempty(strfind(upper(f(j).name),'K27AC'))
            continue;
        end;
        K27{end+1}=[d '/' f(j).name];
    end;
end;

for i=1:length(K27),
    A=load(K27{i});
    for i=1:length(A.S),
	   peak = A.S{i}; 
       bitMap.(peak.chr)(peak.from:peak.to) = bitset(bitMap.(peak.chr)(peak.from:peak.to),4);
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


mean(bitget(bitMap.chr1,1))

% mask H3K4me1 peaks


% bitmap
% 1 = N
% 2 = within gene (up to 15K, including introns)
% 3 = p300 peak (within any condition)
% 4 = H3K27ac peak (within any condition)
% 5 = Unique region (100% unique, >=2 mismatches, read length of 50)

% for example - this could serve as a negative control for enhancer regions
% no N's, not near/in genes, no known enhancer marks, and unique (34% of genome?)
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
    A = ~bitget(bitMap.(c),1) & ~bitget(bitMap.(c),2) & bitget(bitMap.(c),3) & bitget(bitMap.(c),4) & bitget(bitMap.(c),5);
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
