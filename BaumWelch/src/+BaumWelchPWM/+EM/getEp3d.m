% this is used in non all to all mode in backward algorithm, used for the E step after PWM step
% TODO: might be inefficient, many recalculation due to length repetition. could be optimized
% Ret - N x m x length(ts)
function ret = getEp3d(theta, params, X, ts, kronMN, matSize)
    % N x m x length(ts)
    [N, L] = size(X);
    ret = -inf(N, params.m, length(ts));
    for i = 1:length(ts)
        if ts(i) > L
            continue;
        end
        ret(:, :, i) = BaumWelchPWM.EM.getEp(theta, params, X, ts(i), kronMN, matSize);
    end
end