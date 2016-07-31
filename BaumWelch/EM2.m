
function [startT, T, E, likelihood] = EM2(X, m, n, itter)
    % X - 1 x k emission variables
    % m - ammount of possible states (y)ee
    % n - amount of possible emmissions (x)
    % initial estimation parameters
    % [startT, T, E] = genRandParams(m, n);
    % [T, E] = hmmtrain(X,T,E);
    % T(:, 1) = [];
    % startT = T(1, :);
    % T(1, :) = [];
    % E(1, :) = [];
    % return;
    
    k = length(X);
    epsilon = 10 ^ -4;
    [startT, T, E] = genRandParamsExt2(m, n);
    % errors = [];
    likelihood = [];
    itters = 1:itter;
    for i=itters;

        % m x k
        [alpha, scale] = forwardAlg2nd(X, startT, T, E);
        % m x k
        beta = backwardAlg2nd(X, startT, T, E, scale);
        % m x k
        % gamma_t(i) = P(y_t = i|x_1:k)
        gamma = alpha .* beta ./ repmat(sum(alpha .* beta, 1), [m, 1]);
        
        xi = zeros(m);
        for t = 1 : k - 1
            newXi = (alpha(:, t) * (beta(:, t + 1) .* E(:, X(t), X(t + 1)))') .* T;
            xi = xi + newXi / sum(sum(newXi));
        end

        % update estimation parameters
        startT = gamma(:, 1);
        T =  xi;
        for i = 1 : n
            for j = 1 : n
                ijCouple = ([0, X] == i & [X, 0] == j);
                E(:, i, j) = E(:, i, j) + sum(gamma(:, ijCouple(1:end-1)), 2);
            end 
        end
        for i = 1 : m
            for j = 1 : n
                E(i, j, :) = E(i, j, :)  / sum(E(i, j, :), 3);
            end 
        end

        T = bsxfun(@times, T, 1 ./ sum(T, 2));
        likelihood(end + 1) = sum(log(scale));
        if length(likelihood)>1 & abs((likelihood(end) - likelihood(end -1)) / likelihood(end)) < epsilon
            % likelihood converged
            likelihood = likelihood(end);
            return
        end
    end
end




function [startT, T, E] = genRandParamsExt2(m, n)
    % normalized random probabilities
    startT = rand(m, 1);
    startT = startT / sum(startT);
    % epsilon = 0.03;
    % m x m 
    % T = eye(m) + epsilon * rand(m);
    T = rand(m);
    T = bsxfun(@times, T, 1 ./ sum(T, 2));
    % m x n
    E = rand(m, n, n);
    for i = 1 : m
        for j = 1 : n
            E(i, j, :) = E(i, j, :)  / sum(E(i, j, :), 3);
        end 
    end
end


% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n x n emission matrix E_ij means x_t = j | y_t = i
% X - 1 x k emission variables
function [alpha, scale] = forwardAlg2nd(X, startT, T, E)
    % alpha(i,j) P(y_j=i| x_1, ...x_j, startT, T, E)
    % scale(i) = P(x_i| startT, T, E)
    % meaning: alpha_j(i) = alpha(j, i)
    % m x k
    k = length(X);
    [m, n, ~] = size(E);

    alpha = zeros(m, k);
    scale = zeros(1, k);

    % P(x_1=j | y_1=i)
    EStart(:, :) = sum(E, 2);
    EStart = EStart ./ repmat(sum(EStart, 2), [1, n]);
    alpha(:, 1) = startT .* EStart(:, X(1));
    scale(1) = sum(alpha(:, 1));
    for t = 2:k
        % newAlphas = E(:, X(t)) .* sum(repmat(alpha(:, t-1), [1, m]) .* T, 1).';
        newAlphas = (T.' * alpha(:, t-1)) .* E(:, X(t-1), X(t));
        % (T.' * alpha(:, t-1)) is a vector, where v_i is P(y_t=i|x_1...x_t-1)
        scale(t) = sum(newAlphas);
        alpha(:, t) = newAlphas / scale(t);
    end
    % checked
end


% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% X - k x 1 emission variables
function beta = backwardAlg2nd(X, startT, T, E, scale)

    % beta(i,t) P( x_t+1, ...x_k| y_t=i, startT, T, E)
    % beta(t) is the t'th column of beta
    % meaning: beta_j(i) = beta(j, i)

    k = length(X);
    m = length(T);
    % m x k

    beta = ones(m, k);

    for t = k - 1 : -1 : 1
        % newBetas = sum(repmat(E(:, X(t+1)) .* beta(:, 1), [1, m]).' .* T, 2);
        % beta= cat(2, newBetas, beta);
        % beta(:, t) = (T * (E(:, X(t+1), X(t+2)) .* beta(:, t+1))) / scale(t);
        beta(:, t) = (T * (E(:, X(t), X(t+1)) .* beta(:, t+1))) / scale(t);
    end

    % checked
end
