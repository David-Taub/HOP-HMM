% function out = sumDim(M, dims)
%     out = M;
%     for i = 1 : length(dims)
%         out = sum(out, dims(i));
%     end
% end

function out = mulDim(M, dims)
    Msize = size(M);
    nonSummed = 1 : length(Msize);
    nonSummed(dims) = [];
    Mper = permute(M, [nonSummed, dims]);
    Mres = reshape(Mper, [Msize(nonSummed), prod(Msize(dims))]);
    out = prod(Mres, length(nonSummed)+1);
end