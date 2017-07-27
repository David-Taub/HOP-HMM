% this is used in non all to all mode in backward algorithm, used for the E step after PWM step
% TODO: might be inefficient, many recalculation due to length repetition. could be optimized
function ret = getEp3d(theta, params, Xs, ts, kronMN, matSize)
    % N x m x k
    ret = zeros(params.N, params.m, params.k);
    for i = 1:params.k
        ret(:, :, i) = BaumWelchPWM.getEp(theta, params, Xs, ts(i), kronMN, matSize);
    end
end