function [G] = read_gene_information
if ~exist('G','var'),
    if exist('/cs/cbio/tommy/Enhancers/Data/knownGeneXref.mm9.mat','file'),
	load('/cs/cbio/tommy/Enhancers/Data/knownGeneXref.mm9.mat');
    else
	% flyBaseGene.txt
	fid = fopen('/cs/cbio/tommy/Enhancers/Data/knownGeneXref.mm9.txt'); tmp = textscan(fid, '%s%s%c%d%d%d%d%s%s%*[^\n]','Delimiter','\t','bufsize',32767); fclose(fid);
	% acc chr strand txnst txnen cdsst cdsen name desc
	for i=1:length(tmp{1})
	    G(i).acc = tmp{1}{i};
	    G(i).chr = tmp{2}{i};
	    G(i).strand = tmp{3}(i);
	    if G(i).strand  == '+',
		G(i).TSS = tmp{4}(i); G(i).TES = tmp{5}(i); G(i).CDS = tmp{6}(i); G(i).CDE = tmp{7}(i);
	    else
		G(i).TSS = tmp{5}(i); G(i).TES = tmp{4}(i); G(i).CDS = tmp{7}(i); G(i).CDE = tmp{6}(i);
	    end
	    G(i).name = tmp{8}{i};
	    G(i).desc = tmp{9}{i};
	end
	clear tmp;
	save('knownGeneXref.mm9.mat','G');
    end
else
    load('knownGeneXref.mm9.mat');
end
