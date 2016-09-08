function out = sumDim(M, dims)
    out = M;
    for i = 1 : length(dims)
        out = sum(out, dims(i));
    end
end