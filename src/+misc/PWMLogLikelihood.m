

function res = PWMLogLikelihood(params, Xs1H, pwmId)
    res = PWMLogLikelihoodFor(params, Xs1H, pwmId);
    % res = PWMLogLikelihoodConv(params, Xs1H, pwmId)
end

% res - N x L
% Xs1H - N x L x n
function res = PWMLogLikelihoodConv(params, Xs1H, pwmId)
    % PWMs - k x n x J
    pwm = params.PWMs(pwmId, :, :);
    pwm(1, :, params.lengths(pwmId)+1:end) = 1;
    N = size(Xs1H, 1);
    flippedLogPWM = log(permute(pwm(end:-1:1, end:-1:1, end:-1:1), [1, 3, 2]));
    res = convn(Xs1H, flippedLogPWM, 'valid');
    %todo: check this -1 in J-1
    res = [res, -inf(N, params.J-1)];
    res = [-inf(N, params.J-1), res];
end



% res - N x L
% Xs1H - N x L x n
function res = PWMLogLikelihoodFor(params, Xs1H, pwmId)
    % PWMs - k x n x J
    N = size(Xs1H, 1);
    L = size(Xs1H, 2);
    res = ones(N, L) .* eps;
    % 1 x n x length -> 1 x length x n
    pwm = permute(params.PWMs(pwmId, :, 1:params.lengths(pwmId)), [1, 3, 2]);
    % 1 x length x n -> N x length x n
    pwm = repmat(pwm, [N, 1, 1]);
    for i = 1:L - params.lengths(pwmId) + 1
        res(:, i) = prod(sum(pwm .* Xs1H(:, i:i+params.lengths(pwmId) - 1, :), 3), 2);
    end
    res = log(res);
end