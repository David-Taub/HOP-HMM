

function [vals, inds] = maxN(A, d, n)
    vals = [];
    inds = [];
    for i = 1:n
        [new_vals, new_inds] = max(A, [], d);
        A(sub2ind(size(A), [1:size(A, 1)], new_inds')) = -inf;
        vals = cat(d, vals, new_vals);
        inds = cat(d, inds, new_inds);
    end
end

