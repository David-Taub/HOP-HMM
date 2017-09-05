% res - N x k
% PWMsRep - N x J x n x k
% mask - N x J x 1 x k
% Xs1H - N x L x n
function res = getPWMp(J, PWMsRep, Xs1H, t, mask)

    % N x J x n
    lastJXs1H = Xs1H(:, t:t+J-1, :);
    k = size(PWMsRep, 4);
    % N x J x n x k
    res = PWMsRep .* repmat(lastJXs1H, [1, 1, 1, k]);
    % N x J x n x k -> % N x J x 1 x k
    res = sum(res, 3);
    res = res + mask;
    % N x J x 1 x k -> N x 1 x 1 x k
    res = sum(log(res), 2);
    % N x 1 x 1 x k -> N x k
    res = permute(res, [1, 4, 2, 3]);
end