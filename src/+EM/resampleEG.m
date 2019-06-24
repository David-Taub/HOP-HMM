
function [E, G] = resampleEG(params, E, G)
    EFlat = E(:, :);
    randTheta = misc.genTheta(params, false);
    E(:, :) = replaceSimRows(E(:, :), randTheta.E(:, :));
    G(:, :) = replaceSimRows(G(:, :), randTheta.G(:, :));
end


% v - m x 1
function ret = findLowVals(v)
    THRESHOLD = 2.5;
    zScores = v;
    for i = 1:length(v)
        vWithoutI = [v(1:i-1), v(i+1:end)];
        zScores(i) = (mean(vWithoutI, 2) - v(i)) ./ std(vWithoutI);
    end
    ret = zScores > THRESHOLD;
end


% M - m1, m2
function ret = findSimilarRows(M)
    v = pdist(M);
    simMask = findLowVals(v);
    ret = squareform(simMask);
end


% M - m1, m2
function M = replaceSimRows(M, randM)
    for i = size(M, 1) + 1;
        simMask = findSimilarRows(M);
        rowToReplace = find(any(simMask > 0, 2), 1);
        if length(rowToReplace) == 0
            return;
        end
        fprintf('Resampling row %d in matrix %d x %d\n', i, size(M, 1), size(M, 2));
        M(rowToReplace, :) = randM(rowToReplace, :);
    end
    error('Should not reach here, replaced too many rows. Bad randM?');
end