
function [alpha, scale] = forwardAlg(X, startT, T, E)
    if length(size(E)) == 2
        [alpha, scale] = forwardAlg1st(X, startT, T, E);
    else
        [alpha, scale] = forwardAlg2nd(X, startT, T, E);
    end
end


% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% X - 1 x k emission variables
function [alpha, scale] = forwardAlg1st(X, startT, T, E)
    % alpha(i,j) P(y_j=i| x_1, ...x_j, startT, T, E)
    % scale(i) = P(x_i| startT, T, E)
    % meaning: alpha_j(i) = alpha(j, i)
    % m x k
    k = length(X);
    m = length(T);
    alpha = zeros(m, k);
    scale = zeros(1, k);
    alpha(:, 1) = startT .* E(:, X(1));
    scale(1) = sum(alpha(:, 1));
    for t = 2:k
        % newAlphas = E(:, X(t)) .* sum(repmat(alpha(:, t-1), [1, m]) .* T, 1).';
        newAlphas = (T.' * alpha(:, t-1)) .* E(:, X(t));
        scale(t) = sum(newAlphas);
        alpha(:, t) = newAlphas / scale(t);
    end
    % checked
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