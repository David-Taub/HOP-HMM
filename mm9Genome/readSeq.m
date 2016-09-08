% reads .seq format file (Tommy's format) and returns the sequence as numbers
function out = readSeq(filePath, L)
    % matPath = strcat(filePath, '.mat');
    % if exist(matPath, 'file') == 2
    %     out = load(matPath);
    %     out = out.out;
    %     return;
    % end
    fid = fopen(filePath);
    seqsCells = textscan(fid, '%s%s%*[^\n]','delimiter','\t');
    offset = -floor(L/2) + 1: floor(L/2);
    seqs = cellfun(@(x) nt2int(upper(x(ceil(length(x)/2 + offset)))), seqsCells{2}, 'Un', 0);
    out = cell2mat(seqs);
end
