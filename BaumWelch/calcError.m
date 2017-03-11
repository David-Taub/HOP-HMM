
function e = calcError(Y, YEst)
    m = max([Y,YEst]);
    Y = [1:m , Y];
    YEst = [1:m , YEst];
    ct = crosstab(Y, YEst) - eye(m);
    correct = 0;
    matches = [];
    for i = 1:m
        [v1, i1] = max(ct, [], 1);
        [v2, i2] = max(v1, [], 2);
        ct(:, i2) = -inf;
        ct(i1(i2), :) = -inf;
        % match is the mapping between the orignal states numbers and the estimated states number
        % the state number itself has no meaning, only the prob. and the states order is important 
        matches(:, end + 1) = [i1(i2), i2];
        correct = correct + v2;
    end
    e = 1 - (sum(correct) / length(Y));
end