% slice - N x m x k
% Y - m x k transfer probability matrix between mode bases and their PWM modes
% pcPWMp - N x k x L-1 - precomputed likelihood of the sequences and the PWM.
%            at place i, j , t we have the likelihood of seq_i at t+1, with PWM j
% res - N x m
%
function res = PWMstep(slice, Y, ts, pcPWMp, J)
    [N, m, k] = size(slice);
    % N x m x k -> m x k x N
    sliceRep = permute(slice, [2,3,1]);
    % probability to get into the PWM submode
    % N x k x m
    R = permute(bsxfun(@times, sliceRep, Y), [3,2,1]);
    % N x k
    % PWMProb = BaumWelchPWM.getPWMp(PWMsRep, Xs1H, t);
    % in the forward algorithm we ask for t - lengths, and this is a fix for it to be 0
    subscripts = [repmat(1:N, [1, k]); kron([1:k; ts' + J], ones(1, N))];
    indices = matUtils.matSub2ind(size(pcPWMp), subscripts);
    % N x k
    tsPWMp = reshape(pcPWMp(indices), [N, k]);
    % N x k x m, N x k
    res = matUtils.sumDim(bsxfun(@times, R, tsPWMp), 2);
end