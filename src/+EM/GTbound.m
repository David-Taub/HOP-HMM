% this code is very dirty it is built out of heuristics that help the algorithm converge to the correct result
% G - m x k
% T - m x m
function [G, T, startT] = GTbound(params, G, T, doResample)
    if params.m == 1
        return;
    end
    G = exp(G);
    T = exp(T);
    originG = G;
    originT = T;


    params.maxT < T
    params.minT > T
    sum(params.maxG < G, 2)
    sum(params.minG > G, 2)

    for i = 1:5
        T = max(T, params.minT);
        G = max(G, params.minG);
        T = min(T, params.maxT);
        G = min(G, params.maxG);
        s = sum(G, 2) + sum(T, 2);

        T = T ./ repmat(s, [1, params.m]);
        G = G ./ repmat(s, [1, params.k]);
    end

    if doResample
        G = makeDifferent(params, G);
    end
    % [G, T] = balanceGTweights(params, G, T);
    % T = limitTDiag(params, T);
    T = limitCrossEnhancers(params, T, params.PCrossEnhancers);

    % force always start in one of the backgrounds modes
    startT = [ones(params.m - params.backgroundAmount, 1) * eps; ones(params.backgroundAmount, 1) * (1 - (eps * (params.m - params.backgroundAmount))) / params.backgroundAmount];
    fprintf('%d/%d (G) and %d/%d (T) changed in bounding process\n', sum(sum(originG ~= G, 2), 1), params.m * params.k, sum(sum(originT ~= T, 2), 1), params.m ^ 2);
    G = log(G);
    T = log(T);
    startT = log(startT);
end

function G = makeDifferent(params, G)

    [m, k] = size(G);
    for i = 1:m
        for j = i+1:m
            n1 = norm(G(i, :));
            n2 = norm(G(j, :));
            similarity = (G(i, :) ./ n1) * (G(j, :)./ n2)';
            if similarity > 0.95
                fprintf('state %d and %d are too similar (%.2f), resampling %d\n', i, j, similarity, j);
                G(j, :) = rand(1, k);
                G(j, :) = G(j, :) ./ (sum(G(j, :), 2) * params.PTotalBaseToSub);
            end
        end
    end
end

% replaced with min max prob values scheme
% % T - m x m
% function T = limitTDiag(params, T)
%     for i = 1 : params.m - params.backgroundAmount
%         for t = 0 : params.backgroundAmount - 1
%             % fix enhancer length
%             T = makeDiagonalDominant(params, T, i, params.m - t, params.PEnhancerToBackground);
%             % fix background length
%             T = makeDiagonalDominant(params, T, params.m - t, i, params.PBackgroundToEnhancer);
%         end
%         for j = 1 : params.m - params.backgroundAmount
%             % limits cross enhancers occurrences
%             T = makeDiagonalDominant(params, T, i, j, params.PCrossEnhancers);
%         end
%     end
% end

% % T - m x m
% % transfers from T(i,j) to T(i,i)
% function T = makeDiagonalDominant(params, T, i, j, threshold)
%     if isExceedThreshold(params, threshold, T(i, j)) & i ~= j
%         T(i, i) = T(i, i) + T(i, j) - threshold;
%         T(i, j) = threshold;
%     end
% end


% T - m x m
% transfers from T(i,j) to T(i,i)
function T = limitCrossEnhancers(params, T, threshold)
    for i = 1:params.m - params.backgroundAmount
        for j = 1:params.m - params.backgroundAmount
            if j == i
                continue;
            end
            if T(i, j) > threshold
                T(i, params.m) = T(i, params.m) + T(i, j) - threshold;
                T(i, j) = threshold;
            end
        end
    end
end

% % G - m x k
% % T - m x m
% function [G, T] = balanceGTweights(params, G, T)
%     for i = 1:params.m - 1
%         if isExceedThreshold(params, params.PTotalBaseToSub, sum(G(i, :), 2))
%             G(i, :) = G(i, :) .* (params.PTotalBaseToSub / sum(G(i, :), 2));
%             T(i, :) = T(i, :) .* ((1-params.PTotalBaseToSub) / (sum(T(i, :), 2) + eps));
%         end
%     end
%     % make background with eps sub modes (motifs)
%     G(params.m, :) = G(params.m, :) .* ((params.k*eps) / sum(G(params.m, :), 2));
%     T(params.m, :) = T(params.m, :) .* ((1-(params.k*eps)) / (sum(T(params.m, :), 2) + eps));
% end

% function res = isExceedThreshold(params, thresh, val)
%     res = val > thresh * params.maxPRatio || val < tresh / params.maxPRatio;
%     % res = abs(val - thresh) > thresh * params.maxPRatio;
% end

