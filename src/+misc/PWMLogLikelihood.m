

function res = PWMLogLikelihood(PWMs, lengths, Xs1H, pwmId)
    res = PWMLogLikelihoodFor(PWMs, lengths, Xs1H, pwmId);
    % res = PWMLogLikelihoodConv(PWMs, lengths, Xs1H, pwmId)
end

% res - N x L
% Xs1H - N x L x n
function res = PWMLogLikelihoodConv(PWMs, lengths, Xs1H, pwmId)
    % PWMs - k x n x J
    J = size(PWMs, 3);
    pwm = PWMs(pwmId, :, :);
    pwm(1, :, lengths(pwmId)+1:end) = 1;
    N = size(Xs1H, 1);
    flippedLogPWM = log(permute(pwm(end:-1:1, end:-1:1, end:-1:1), [1, 3, 2]));
    res = convn(Xs1H, flippedLogPWM, 'valid');
    %todo: check this -1 in J-1
    res = [res, -inf(N, J-1)];
    res = [-inf(N, J-1), res];
end



% res - N x L
% Xs1H - N x L x n
% PWMs - k x n x J
function res = PWMLogLikelihoodFor(PWMs, lengths, Xs1H, pwmId)
    N = size(Xs1H, 1);
    L = size(Xs1H, 2);
    res = ones(N, L) .* eps;
    % 1 x n x length -> 1 x length x n
    pwm = permute(PWMs(pwmId, :, 1:lengths(pwmId)), [1, 3, 2]);
    % 1 x length x n -> N x length x n
    pwm = repmat(pwm, [N, 1, 1]);
    for i = 1:L - lengths(pwmId) + 1
        res(:, i) = prod(sum(pwm .* Xs1H(:, i:i + lengths(pwmId) - 1, :), 3), 2);
    end
    res = log(res);
end