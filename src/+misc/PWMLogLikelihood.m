
% res - N x L
% Xs1H - N x L x n
function res = PWMLogLikelihood(params, Xs1H, pwmId)
    % PWMs - k x n x J
    pwm = params.PWMs(pwmId, :, :);
    pwm(1, :, params.lengths(pwmId)+1:end) = 1;
    N = size(Xs1H, 1);
    flippedLogPWM = log(permute(pwm(end:-1:1, end:-1:1, end:-1:1), [1, 3, 2]));
    res = convn(Xs1H, flippedLogPWM, 'valid');
    %todo: check this -1 in J-1
    res = [res, -inf(N, params.J-1)];
end