
% res - N x L
% onePaddedPWM - 1 x n x J
% Xs1H - N x L x n
function res = PWMLogLikelihood(params, Xs1H, pwmId)
    % PWMs - k x n x J
    pwm = params.PWMs(pwmId, :, :);
    pwm(1, :, params.lengths(pwmId)+1:end) = 1;
    N = size(Xs1H, 1);
    flippedLogPWM = log(permute(pwm(end:-1:1, end:-1:1, end:-1:1), [1, 3, 2]));
    res = convn(Xs1H, flippedLogPWM, 'valid');
    res = [res, -inf(N, params.J)];
end