function ret = logMatSum(A, dim)
    ret = log(sum(exp(A), dim));
end
