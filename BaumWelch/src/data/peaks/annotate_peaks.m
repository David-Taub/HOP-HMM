function [S] = annotate_peaks(S)
N=length(S);

% read gene information
global G; if isempty(G), G=read_gene_information; end;
Gchr = {G.chr}; Gnames = {G.name};

Ybnd = [S{:}];
Pchrs = {Ybnd.chr};
chrs = unique(Pchrs);

% add center point
for i=1:length(S), p=S{i}; c=p.chr; f=p.from; t=p.to; S{i}.pos=ceil(mean(f,t)); end

% go over chromosomes
% chrs = {'chr3LHet'});
for ic=1:length(chrs),
    chr = chrs{ic},;
    % find genes
    GI = find(ismember(Gchr,chr));
    if isempty(GI), continue; end
    Gstart = [G(GI).TSS]'; Gend = [G(GI).TES]'; GGpos = [Gstart,Gend]; Gstrand = [G(GI).strand]=='+';

    % find peaks
    II = find(ismember(Pchrs,chr));
    % go over peaks
    for i=1:length(II),
	pid = II(i); p=S{pid};

	% find nearest TSS (just one)
	[Gdist,j]=min(abs(p.pos-Gstart)); gid=GI(j); Gdist=p.pos-Gstart(j); Gname=G(gid).name; Gdesc=G(gid).desc;
	if G(gid).strand=='-', Gdist=-Gdist; end;
	p.gTSS = Gname; p.gTSSdesc = Gdesc; p.gTSSdist = Gdist; p.gTSSid=gid;

	% find nearest TES (just one)
	[Gdist2,j]=min(abs(p.pos-Gend)); gid=GI(j); Gdist2=p.pos-Gend(j);  Gname=G(gid).name; Gdesc=G(gid).desc;
	if G(gid).strand=='-', Gdist2=-Gdist2; end;
	p.gTES = Gname; p.gTESdesc = Gdesc; p.gTESdist = Gdist2; p.gTESid=gid;

	% find 3 nearest genes 
	dists = minabs(p.pos-GGpos); [~,J]=sort(abs(dists));
	GG = {G(GI(J)).name}; Gacc = {G(GI(J)).acc};
	Gtmp = GG; Gtmpd = dists(J)'; Gi = {G(GI(J)).acc}; 
	Genes = cell(1,3); Gdists = NaN * ones(1,3); GeneAcc = cell(1,3); Gids = NaN * ones(1,3);
	for r=1:3,
	    if ~isempty(Gtmp),
		ii=GI(J(1));
		Genes(r) = Gtmp(1); Gdists(r) = Gtmpd(1); GeneAcc{r} = Gi{1}; Gids(r)=ii;
		if G(ii).strand=='-', Gdists(r) = -Gdists(r); end;
		% remove all other copies of the same gene name
		III = ismember(Gtmp,Gtmp(1));
		Gtmp(III)=''; Gtmpd(III)=[]; J(III)=[]; Gi(III)=[];
	    end
	end
	p.Genes = Genes; p.Gdists = Gdists; p.GeneAcc=GeneAcc; p.GeneIDs = Gids;
	clear Gtmp Gtmpd r III gid GeneAcc;

	% classifications
	if abs(p.gTSSdist)<=500, p.class = 'Promoter';
	elseif p.gTESdist>-1000 & p.gTESdist<3000, p.class = 'Terminal';
	elseif abs(p.gTSSdist)<=5000, p.class = 'Proximal';
	elseif abs(p.gTSSdist)<=50000, p.class = 'Distal';
	else p.class = 'Intergenic'; end

	S{pid} = p;
    end
end
