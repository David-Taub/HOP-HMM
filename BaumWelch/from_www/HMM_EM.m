function [model,log_like] = HMM_EM(O_all,N)
% learn HMM parameters using EM
%
% input:
%   O_all:  1 x seqNum cell, each is a
%               1 x T, observed sequence, with M symbols (1,..,M)
%   N:      number of hidden states
%
% output:
%   model:      a model, which contains the following estimated parameters
%       A:T          N x N, transition matrix, a_ij = Prb(q_j|q_i)
%       B:E          N x M, emission matrix, b_ij = Prb(o_j|q_i)
%       P:          N x 1, prior probabilities
%   log_like:   log likelihood of each iteration

M = max(O_all);     % assume all symbols have occured at least once
conv_prec = 1e-6;
max_iter = 1000;

% random initialization
P = rand(N,1)+eps; P = P/sum(P);
A = rand(N,N)+eps; A = bsxfun(@times,A,1./sum(A,2));
B = rand(N,M)+eps; B = bsxfun(@times,B,1./sum(B,2));

log_like = zeros(max_iter,1);
% errors = [];
for i = 1:max_iter
    new_P = zeros(size(P));
    new_A = zeros(size(A));
    new_B = zeros(size(B));
    O = O_all;
    T = length(O);
    
    % compute forward and backward probabilities
    [alpha,scale_alpha] = compForwardProb(O,A,B,P);
    [beta] = compBackwardProb(O,A,B,scale_alpha);

    % compute posterior probabilities (E-step)
    [gamma] = compGamma(alpha,beta);

    % compute averaged joint posterior (q_i,q_j|O)
    ksi = zeros(N);
    for t = 1:T-1
        ksi_tmp = (alpha(:,t) * (beta(:,t+1).*B(:,O(t+1)))') .* A;
        ksi = ksi + ksi_tmp / sum(sum(ksi_tmp));
    end
    
    % update parameters (M-step)
    new_P = new_P + gamma(:,1);
    new_A = new_A + ksi;
    for k = 1:M
        new_B(:,k) = new_B(:,k) + sum(gamma(:,O==k),2);
    end
    
    % evaluate log-likelihood
    log_like(i) = log_like(i) - sum(log(scale_alpha));
    
    
    % normalize update
    P = new_P;
    A = bsxfun(@times,new_A,1./sum(new_A,2));
    B = bsxfun(@times,new_B,1./sum(new_B,2));
    
    % errors(end + 1) = norm(B(:) - BOrig(:));
    % determine if converged
    if i > 2
        log_like_change = abs(1-log_like(i-1)/log_like(i));
        if log_like_change < conv_prec
            break;      % converged
        end
    end
end
% figure
% plot(errors)

model.A = A;
model.B = B;
model.P = P;
log_like = log_like(1:i);