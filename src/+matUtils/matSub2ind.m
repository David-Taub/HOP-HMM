	
% very similar to sub2ind, but receives the subscripts as matrix
% matSize - 1 x k
% subscripts - k x n
% indices - 1 x n
function indices = matSub2ind(matSize, subscripts)
    subtractedSub = subscripts;
    cumMatSize = cumprod([1, matSize]);
    cumMatSize = cumMatSize(1, 1:end - 1);

    % note: the following two lines is the same as commented line, but much faster
    % subtractedSub(2:end, :) = subtractedSub(2:end, :) - 1;
    subtractedSub = subtractedSub - 1;
    subtractedSub(1,:) = subtractedSub(1,:) + 1;
    
    indices = double(cumMatSize) * double(subtractedSub);
end
