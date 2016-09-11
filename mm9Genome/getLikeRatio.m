
% X - N x L
function likeRatio = getLikeRatio(E, X)
    s = size(E);
    posE = reshape(E(1,:), s(2:end));
    negE = reshape(E(2,:), s(2:end));
    likePos = getLogLikes(posE, X);
    likeNeg = getLogLikes(negE, X);
    likeRatio = likePos ./ likeNeg; %high / low = high
end
    
    
function logLikes = getLogLikes(E, seqs)
    [N, L] = size(seqs);
    order = matDim(E);
    indices = getIndeices1D(seqs, order);
    indices = reshape(indices, [L - order + 1, N]);
    logLikes = sum(log(E(indices)), 1);
end

