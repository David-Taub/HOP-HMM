% PWMsRep - J x L-J+1  x n x k
% X1H - 1 x L x n
% res - L-J+1 x k
% mask - J x L-J+1 x 1 x k
function res = getPWMpMax(J, PWMsRep, X1H, t, mask)
    k = size(PWMsRep, 4);
    [~, L, n] = size(X1H);

    % J x L-J+1 x n
    slidingWinX1H = zeros(J, L - J + 1, n);
    for i = 1:J
        slidingWinX1H(i, :, :) = X1H(1, i:end - J + i, :);
    end
    % J x L-J+1 x n x k
    slidingWinX1H = repmat(slidingWinX1H, [1,1,1,k]);
    res = PWMsRep .* slidingWinX1H;
    res = sum(res, 3);
    res = res + mask;
    res = prod(res, 1);
    res = permute(res, [2,4,1,3]);
end