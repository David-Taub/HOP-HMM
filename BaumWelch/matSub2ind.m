
% very similar to sub2ind, but receives the subscripts as matrix
% matSize - 1 x k
% subscripts - k x n
% indices - 1 x n
function indices = matSub2ind(matSize, subscripts)
    subtractedSub = subscripts.';
    subtractedSub(2:end, :) = subtractedSub(2:end, :) - 1;
    % subtractedSub
    cumMatSize = cumprod([1, matSize]);
    cumMatSize = cumMatSize(1, 1:end - 1);
    indices = cumMatSize * subtractedSub;
end
