% this is used in non all to all mode in backward algorithm, used for the E step after PWM step
% TODO: might be inefficient, many recalculation due to length repetition. could be optimized
% Ret - N x m x length(ts)
function ret = getEp3d(theta, params, Xs, ts, kronMN, matSize)
    % N x m x length(ts)
    ret = zeros(params.N, params.m, length(ts));
    for i = 1:length(ts)
        ret(:, :, i) = BaumWelchPWM.getEp(theta, params, Xs, ts(i), kronMN, matSize);
    end
end