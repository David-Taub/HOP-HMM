% this is used in non all to all mode in backward algorithm, used for the E step after PWM step
% TODO: might be inefficient, many recalculation due to length repetition. could be optimized
% Ret - N x m x length(ts)
function ret = getEp3d(theta, params, X, ts)
    % N x m x length(ts)
    [N, L] = size(X);
    matSize = [params.m , params.n * ones(1, params.order)];
    kronMN = kron(1:params.m, ones(1, N));
    ret = -inf(N, params.m, length(ts));
    for i = 1:length(ts)
        if ts(i) > L
            continue;
        end
        ret(:, :, i) = EM.getEp(theta, params, X, ts(i), kronMN, matSize);
    end
end