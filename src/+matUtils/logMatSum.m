function ret = logMatSum(A, dim)
    assert(size(A, dim) > 0);
    % ret = logMatSum2(A, dim, 1, size(A, dim));
    % ret = logMatSum1(A, dim);
    ret = logMatSum3(A, dim);
end

% faster, but calculation error should be bigger
function ret = logMatSum1(A, dim)
    sizeA = size(A);
    sizeRet = sizeA;
    sizeRet(dim) = 1;
    ret = -inf(sizeRet);
    for i = 1:size(A, dim)
        ret = matUtils.logAdd(ret, arraySlice(A, i, dim));
    end
end

% faster, but calculation error should be bigger
function ret = logMatSum3(A, dim)
    Amax = max(A, [], dim);
    repSize = ones(1, length(size(A)));
    repSize(dim) = size(A, dim);
    Asub = A - repmat(Amax, repSize);
    ret = Amax + log(sum(exp(Asub), dim));
    ret(Amax == -inf) = -inf;
    assert(not(any(isnan(ret(:)))));
end

% slower, but calculation error should be smaller
function ret = logMatSum2(A, dim, startIndex, endIndex)
    if startIndex == endIndex
        ret = arraySlice(A, startIndex, dim);
        return
    end
    diff = endIndex - startIndex;
    right = logMatSum2(A, dim, startIndex, startIndex + floor(diff/2));
    left = logMatSum2(A, dim, startIndex + floor(diff/2) + 1, endIndex);
    ret = matUtils.logAdd(right, left);
end

function S = arraySlice(A, I, d)
    if length(size(A)) == 3
        if d == 1
            S = A(I, :, :);
        elseif d == 2
            S = A(:, I, :);
        else
            S = A(:, :, I);
        end
        return;
    end
    s1 = size(A); s2 = size(A); s3 = size(A);
    s1(d) = I-1;
    s2(d) = 1;
    s3(d) = s3(d)-I;
    S = zeros(s2);
    S(:) = A(cat(d, false(s1), true(s2), false(s3)));
end
