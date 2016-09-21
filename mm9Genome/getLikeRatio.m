
% X - N x L
function likeRatio = getLikeRatio(E, X)
    s = size(E);
    if length(s) > 2
        posE = reshape(E(1,:), s(2:end));
        negE = reshape(E(2,:), s(2:end));
    else
        posE = E(1,:);
        negE = E(2,:);
    end
    likePos = getLogLikes(posE, X);
    likeNeg = getLogLikes(negE, X);
    likeRatio = likePos - likeNeg; %high / low = high
end
    
