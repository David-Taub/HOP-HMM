% res - N x k
% PWMsRep - N x J x n x k
function res = getPWMp(J, PWMsRep, Xs1H, t, mask)
    % N x J x n
    lastJXs1H = Xs1H(:,t:t+J-1,:);
    k = size(PWMsRep, 4);
    res = PWMsRep .* repmat(lastJXs1H, [1,1,1,k]);
    res = sum(res, 3);
    res = res + mask;
    res = prod(res, 2);
    res = permute(res, [1,4,2,3]);
end