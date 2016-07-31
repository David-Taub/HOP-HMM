function beta = backwardAlg(X, startT, T, E, scale)
    if length(size(E)) == 3
        beta = backwardAlg2nd(X, startT, T, E, scale);
    else
        beta = backwardAlg1st(X, startT, T, E, scale);
    end
end 

% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% X - k x 1 emission variables
function beta = backwardAlg1st(X, startT, T, E, scale)

    % beta(i,j) P( x_j+1, ...x_k| y_j=i, startT, T, E)
    % beta(t) is the t'th column of beta
    % meaning: beta_j(i) = beta(j, i)

    k = length(X);
    m = length(T);
    % m x k

    beta = ones(m, k);

    for t = k - 1 : -1 : 1
        % newBetas = sum(repmat(E(:, X(t+1)) .* beta(:, 1), [1, m]).' .* T, 2);
        % beta= cat(2, newBetas, beta);
        beta(:, t) = (T * (E(:, X(t+1)) .* beta(:, t+1))) / scale(t);
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