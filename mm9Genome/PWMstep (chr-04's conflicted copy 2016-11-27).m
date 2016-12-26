% slice - N x m * k from alpha or beta of the forward or backwards algorithm 
% PWMsRep - N x J x n x k: N replication of emission matrix of m Jaspar PWMs with length J 
%        true length of i'th PWM< J and given in lengths(i) if a PWM is 
%        shorter than j, it is aligned to the end of the 3rd dimension.
% T - m x m+k transfer probability matrix between mode bases to other modebases and sub modes
% Xs1H - N x J+L x n: one hot map of the sequences, padded with a J zeros prefix
% res - N x m
function res = PWMstep(slice, PWMsRep, Xs1H, T, t)
    J = size(PWMsRep, 2);
    R = bsxfun(@times, slice, shiftdim(Y, -1));
    % N x J x n
    lastJXs1H = Xs1H(:, t+1:t+J, :);
    % N x J x n x k 
    P = bsxfun(@times, PWMsRep, lastJXs1H);
    % N x k
    res = bsxfun(@times, R, sumDim(P, [2, 3]));
end
