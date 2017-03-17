
% remove short seqs
% ACGT -> 1234
function seqs = regularSeqs(seqs, L)
    lengths = cell2mat(cellfun(@size, seqs, 'uni', false));
    centers = repmat(ceil(lengths(:, 2) ./ 2, [1, L]));
    offsets = repmat([-floor(L/2) + 1: floor(L/2)], [length(lengths), 1]);
    keyboard

    seqs = seqs(:, centers + offsets);
    seqs = nt2int(upper(seqs));
    % in some of the data we have NNN which can be any nucleotide
    % out(out==15) = 1;
end
