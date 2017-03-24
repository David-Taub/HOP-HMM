% slice - N x m x k
% Y - m x k transfer probability matrix between mode bases and their PWM modes
% pcPWMp - N x k x L-1 - precomputed likelihood of the sequences and the PWM. 
%            at place i, j , t we have the likelihood of seq_i at t+1, with PWM j
% res - N x m
% 
function res = PWMstep(slice, Y, ts, pcPWMp)
    % N x m x k -> m x k x N
    [N, m, k] = size(slice);
    sliceRep = permute(slice, [2,3,1]);
    % probability to get into the PWM submode
    % N x k x m
    R = permute(bsxfun(@times, sliceRep, Y), [3,2,1]);
    % N x k
    % PWMProb = BaumWelchPWM.getPWMp(PWMsRep, Xs1H, t);
    % N x k x m, N x k
    % in the forward algorithm we ask for t - lengths, and this is a fix for it to be 0
    ts1 = max(ts, 1);
    indices = matUtils.matSub2ind(size(pcPWMp), [repmat(1:N, [1, k]); kron([1:k; ts1'], ones(1, N))]);
    res = matUtils.sumDim(bsxfun(@times, R, reshape(pcPWMp(indices), [N, k])), 2);
    res(ts < 1) = 0;
end