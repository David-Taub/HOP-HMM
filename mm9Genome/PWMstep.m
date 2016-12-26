% slice - N x m * k
% res - N x m
function res = PWMstep(slice, PWMsRep, Xs1H, Y, J, t)

    R = bsxfun(@times, slice, shiftdim(Y, -1));
    % N x J x n
    lastJXs1H = Xs1H(:, t+1:t+J, :);
    % N x J x n x k 
    P = bsxfun(@times, PWMsRep, lastJXs1H);
    % N x k
    res = bsxfun(@times, R, sumDim(P, [2, 3]));
end
