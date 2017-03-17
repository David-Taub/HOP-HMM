
function dim = matDim(M)
    s = size(M);
    s(s==1) = [];
    dim = length(s);
end