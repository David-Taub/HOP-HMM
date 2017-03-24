% reads .seq format file (Tommy's format) and returns the sequence as numbers
function out = readSeq(filePath, L)
    matPath = strcat(filePath, '.mat');
    [~, seqsCells] = textread(filePath, '%s%s%*[^\n]','delimiter','\t','bufsize',20000);
    out = regularSeqs(seqsCells, L);
end

% remove short seqs
% ACGT -> 1234
function out = regularSeqs(seqsCells, L)
    out = zeros(length(seqsCells), L); 
    for i=1:length(seqsCells)
        seqsCells{i}=upper(seqsCells{i}); 
        if length(seqsCells{i}) < L
            continue;
        end
        seqsCells{i}=seqsCells{i}(ceil(length(seqsCells{i})/2) + [-floor(L/2) + 1: floor(L/2)]);
        out(i, :) = nt2int(seqsCells{i});
    end
    out( ~any(out,2), : ) = [];  %remove zero rows
    % in some of the data we have NNN which can be any nucleotide
    out(out==15) = 1;
end
