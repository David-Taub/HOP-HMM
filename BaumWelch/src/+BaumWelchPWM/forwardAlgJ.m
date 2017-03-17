% This model assumes m modes bases that emits like high order HHM,
% where each base can transfer into k submodes that emits with PWM from JASPAR project
% m - sum of enhancer and background modes (not a parameter)
% T - m x m transfer probability matrix between mode bases
% Y - m x k transfer probability matrix between mode bases and their PWM modes
% n - number of alphabet (4, i.e. ACGT)
% startT - m x 1 probabilities of first states
% E - m x n x n x ... x n (order times) alphabet emission probability matrix
% PWMsRep - N x J x n x k: N replication of emission matrix of m Jaspar PWMs with length J 
%        true length of i'th PWM< J and given in lengths(i) if a PWM is 
%        shorter than j, it is aligned to the end of the 3rd dimension.
% lengths - m x 1 length of each motif in the PWM matrix. J = max(lengths)
% Xs - N x L emission variables

function [alpha, scale] = forwardAlgJ(Xs, startT, T, Y, E, PWMsRep, lengths)
    % alpha(N, i, j) P(y_s_j=i| x_s_1, ...x_s_j, startT, T, PWMs)
    % scale(N, i) = P(x_s_i| startT, T, PWMs)
    % m x L
    [N, L] = size(Xs);
    [~, J, n, ~] = size(PWMsRep);
    m = size(T, 1);
    kronMN = kron(1:m, ones(1, N));
    alpha = zeros(N, m, L + J);
    scale = zeros(N, L);
    order = matDim(E) - 1;
    startE = sumDim(E, 2 : order);
    Xs1H = cat(2, zeros(N, J, n), mat23Dmat(Xs, n));
    alpha(:, :, J+1) = (repmat(startT, [1, N]) .* startE(:, Xs(:, 1))).';
    scale(:, 1) = sum(alpha(:, :, J+1), 2);
    matSize = [m , n * ones(1, order)];

    for t = 2:L
        % N x m 
        Ep = getEp(E, Xs, t, m, kronMN, matSize, N, order);
        % N x m 
        newAlphas = (alpha(:, :, t + J - 1) * T) .* Ep;
        % N x m x k
        alphaSlice = alpha(:, :, t + J - lengths - 1);
        newAlphas = newAlphas + PWMstep(alphaSlice, PWMsRep, Xs1H, Y, t);
        % N x 1
        scale(:, t) = sum(newAlphas, 2);
        alpha(:, :, J + t) = bsxfun(@times, newAlphas, 1 ./ scale(:, t));
    end
    alpha = alpha(:, :, J+1:end);
end
