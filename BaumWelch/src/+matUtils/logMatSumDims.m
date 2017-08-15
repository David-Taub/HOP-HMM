function ret = logMatSumDims(A, dims)
    ret = log(matUtils.sumDim(exp(A), dims));
end
