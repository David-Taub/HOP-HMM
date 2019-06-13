% TODO:this code is very dirty it is built out of heuristics
% help the algorithm converge to the correct result
% G - m x k
% T - m x m
function [G, T] = GTbound(params, G, T)
    if params.m == 1
        return;
    end
    G = exp(G);
    T = exp(T);
    originG = G;
    originT = T;
    T(T > params.maxT) = params.maxT(T > params.maxT);
    G(G > params.maxG) = params.maxG(G > params.maxG);
    T(eye(params.m) == 1) = T(eye(params.m) == 1) + 1 - sum(T, 2) - sum(G, 2);
    fprintf('%d/%d (G) and %d/%d (T) changed in bounding process\n', sum(sum(originG ~= G, 2), 1), params.m * params.k, sum(sum(originT ~= T, 2), 1), params.m ^ 2);
    G = log(G);
    T = log(T);
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

