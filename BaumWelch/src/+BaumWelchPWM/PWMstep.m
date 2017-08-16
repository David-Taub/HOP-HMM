% slice - N x m x k
% Gs - m x k transfer probability matrix between mode bases and their PWM modes
% ts - k x 1 indices from which the motifs starts in all N seqs and m states
% pcPWMp - N x k x L-1 - precomputed likelihood of the sequences and the PWM.
%            at place i, j, t we have the likelihood of seq_i at t:t+J, with PWM j
% Eps - N x m x k
% res - N x m x k
function ret = PWMstep(slice, Gs, ts, pcPWMp, Eps, Fs)
    [N, m, k] = size(slice);
    % probability to get into the PWM submode
    % N x m x k
    % in the forward algorithm we ask for t - lengths, and this is a fix for it to be 0
    subZeroMask = ts<1;
    ts(subZeroMask) = 1;
    subscripts = [repmat(1:N, [1, k]); kron([1:k; ts'], ones(1, N))];
    % subscripts has all Nxk indices in [1:N, 1:k, ts(1:k)]
    indices = matUtils.matSub2ind(size(pcPWMp), subscripts);
    % N x 1 x k
    tsPWMp = reshape(pcPWMp(indices), [N, 1, k]);

    % N x 1 x k -> N x m x k
    tsPWMps = repmat(tsPWMp, [1, m, 1]);

    ret = slice + Gs;
    ret = ret + Fs;
    ret = ret + tsPWMps;
    ret = ret + Eps;
    ret = matUtils.logMatSum(ret, 3);
    ret(:,:,subZeroMask) = -inf;
end