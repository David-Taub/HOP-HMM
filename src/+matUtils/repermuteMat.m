% perm - m x 1
function perm = repermuteMat(matOrig, matEst)
    % to vec
    n = size(matOrig, 1);
    % distMat = vectorizedOrig * vectorizedEst';
    distMat = pdist2(matOrig, matEst);
    perm = zeros(n, 1);
    for i = 1:n
        [~, I] = min(distMat(:));
        [I_row, I_col] = ind2sub(size(distMat), I);
        perm(I_row) = I_col;
        distMat(I_row, :) = inf;
        distMat(:, I_col) = inf;
    end
end

