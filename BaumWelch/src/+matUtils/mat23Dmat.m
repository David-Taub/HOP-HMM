
% M - k1 x k2
% M3d - k1 x k2 x n
% example:
% the input M=[1 2 2; 1 3 4], n=4
% will result:
% 1 0 0  0 1 1  0 0 0  0 0 0
% 1 0 0  0 0 0  0 1 0  0 0 1
function M3D = mat23Dmat(M, n)
    M2D = matUtils.vec2mat(M(:).', n);
    M3D = permute(reshape(M2D, [n, size(M)]), [2,3,1]);
    % M3D = repmat(M, [1,1,n]) == repmat(permute(1:n, [3,1,2]), [size(M),1]);
end