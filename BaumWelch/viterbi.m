function Y = viterbi(X, startT, T, E)
    if length(size(E)) == 3
        Y = viterbi2nd(X, startT, T, E);
    else
        Y = viterbi1st(X, startT, T, E);
    end
end
% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% X - k x 1 emission variables
% Y - k x 1 hidden markov variables
function Y = viterbi1st(X, startT, T, E)
    % Y = ViterbiDecode(X,T,E,startT);
    % return;
    % [T, E] = addPrior(startT, T, E);
    % Y = hmmviterbi(X,T,E) - 1;
    % return;
    
    k = length(X);
    m = length(T);
    % m x k
    V = zeros(m,k);
    ind = zeros(m,k);
    % k, 1
    Y = zeros(1, k);

    V(:, 1) = log(startT) + log(E(:, X(1)));
    for t = 2 : k
        [V(:, t), ind(:, t)] = max(repmat(V(:, t-1).', [m, 1]) + log(T.'), [], 2);
        V(:, t) = V(:, t) + log(E(:, X(t)));
    end
    % Y_k = argmax(y, V_k,y)
    [~, Y(k)] = max(V(:, k));
    for t = k-1 : -1 : 1
        Y(t) = ind(Y(t+1), t+1);
    end
    % checked
end


% T - m x m transfer matrix T_ij means y_t = j | y_t-1 = i
% startT - m x 1 probabilities of first states
% E - m x n emission matrix E_ij means x_t = j | y_t = i
% X - k x 1 emission variables
% Y - k x 1 hidden markov variables
function Y = viterbi2nd(X, startT, T, E)
    % Y = ViterbiDecode(X,T,E,startT);
    % return;
    % [T, E] = addPrior(startT, T, E);
    % Y = hmmviterbi(X,T,E) - 1;
    % return;
    
    k = length(X);
    [m, n, ~] = size(E);
    % m x k
    V = zeros(m,k);
    ind = zeros(m,k);
    % k, 1
    Y = zeros(1, k);

    EStart(:, :) = sum(E, 2);
    EStart = EStart ./ repmat(sum(EStart, 2), [1, n]);

    V(:, 1) = log(startT) + log(EStart(:, X(1)));
    for t = 2 : k
        [V(:, t), ind(:, t)] = max(repmat(V(:, t-1).', [m, 1]) + log(T.'), [], 2);
		V(:, t) = V(:, t) + log(E(:, X(t-1), X(t)));
    end
    % Y_k = argmax(y, V_k,y)
    [~, Y(k)] = max(V(:, k));
    for t = k-1 : -1 : 1
        Y(t) = ind(Y(t+1), t+1);
    end
    % checked
end