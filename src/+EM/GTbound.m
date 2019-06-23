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
    T(T < params.minT) = params.maxT(T < params.minT);
    G(G > params.maxG) = params.maxG(G > params.maxG);
    G(G < params.minG) = params.maxG(G < params.minG);
    T(eye(params.m) == 1) = T(eye(params.m) == 1) + 1 - sum(T, 2) - sum(G, 2);
    fprintf('%d/%d (G) and %d/%d (T) changed in bounding process\n', sum(sum(originG ~= G, 2), 1), params.m * params.k, sum(sum(originT ~= T, 2), 1), params.m ^ 2);
    G = log(G);
    T = log(T);
end

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

