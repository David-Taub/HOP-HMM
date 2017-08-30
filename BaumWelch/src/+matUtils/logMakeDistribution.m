function ret = logMakeDistribution(A)
    ret = bsxfun(@minus, A, matUtils.logMatSum(A, length(size(A))));
end
