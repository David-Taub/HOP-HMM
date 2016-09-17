
% very similar to sub2ind, but receives the subscripts as matrix
% matSize - 1 x k
% subscripts - k x n
% indices - 1 x n
function indices = matSub2ind(matSize, subscripts)
    cumMatSize = cumprod([1, matSize]);
    cumMatSize = cumMatSize(1, 1:end - 1);
    indices = double(cumMatSize) * double(subscripts);
    indices = indices - sum(cumMatSize(2:end), 2);
end
