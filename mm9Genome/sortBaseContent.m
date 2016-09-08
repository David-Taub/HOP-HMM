
function sorted = sortBaseContent(seqs)
    tic
    enhCont = [0.23854, 0.2626, 0.25988, 0.23898];
    [N, L] = size(seqs);
    z = zeros(N, 4);
    for i = 1:4
        z(:, i) = sum(seqs == i, 2);
    end
    z = z ./ L;
    score = sum(abs(bsxfun(@minus, z, enhCont)), 2);
    
    [~, scoreInd] = sort(score);
    sorted = seqs(scoreInd, :);
    toc
end