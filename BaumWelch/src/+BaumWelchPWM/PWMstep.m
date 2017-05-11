% slice - N x m x k
% Y - m x k transfer probability matrix between mode bases and their PWM modes
% pcPWMp - N x k x L-1 - precomputed likelihood of the sequences and the PWM.
%            at place i, j , t we have the likelihood of seq_i at t+1, with PWM j
% res - N x m x k
%
function res = PWMstep(slice, Ms, ts, pcPWMp, J)
    [N, m, k] = size(slice);
    % probability to get into the PWM submode
    % N x m x k
    slice = slice .* Ms;
    % in the forward algorithm we ask for t - lengths, and this is a fix for it to be 0
    subscripts = [repmat(1:N, [1, k]); kron([1:k; ts' + J], ones(1, N))];
    indices = matUtils.matSub2ind(size(pcPWMp), subscripts);
    % N x 1 x k
    tsPWMp = reshape(pcPWMp(indices), [N, 1, k]);

    % N x 1 x k -> N x m x k
    tsPWMps = repmat(tsPWMp, [1, m, 1]);
    slice = slice .* tsPWMps;
    % res = sum(slice, 3);
end