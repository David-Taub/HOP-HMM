function ret = logMatSum(A, dim)
    assert(size(A, dim) > 0);
    Amax = max(A, [], dim);
    repSize = ones(1, length(size(A)));
    repSize(dim) = size(A, dim);
    Asub = A - repmat(Amax, repSize);
    ret = Amax + log(sum(exp(Asub), dim));
    ret(Amax == -inf) = -inf;
    assert(not(any(isnan(ret(:)))));
end
