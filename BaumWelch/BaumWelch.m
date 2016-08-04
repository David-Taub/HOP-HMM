% make an EM version that uses context
% using real data of enhancer and non enhancer statistics


% Steps:
% get E T context matrix
% generate text
% EM guess the best E T 

% first version: (didn't really do this, although Tommy said I should)
% get 1st level E T
% viterbi - get states Y
% get couples statistics, build 3D E T

% second version: 
% really change EM file
function BaumWelch()
    close all;
    repeat = 5;
    m = 2;
    n = 4;
    epsilonT = 0.005;
    e1 = []; e2 = [];
    like1 = []; like2 = []; %e3 = [];
    [startTReal, TReal, EReal] = genGenomeParams(epsilonT);
    % [startTReal, TReal, EReal] = genStrong2dParams(epsilonT);

    for i =1:repeat

        k = 10000;
        % Create X and Y
        fprintf('Genarated Params\n');
        % [startTReal, TReal, EReal]
        fprintf('Creating X and Y length: %d\n', k);
        [X, Y] = genHmm2nd(k, startTReal, TReal, EReal);
        % 1
        fprintf('Running EM\n');
        [startTMine, TMine, EMine, likelihood, gamma] = EM(X, m, n, 400);
        like1(end+1) = sum(log(gamma(vec2mat(Y,2))))
        fprintf('Viterbi\n');
        YEst = viterbi(X, startTMine, TMine, EMine);
        err = calcError(Y, YEst);
        e1(end + 1) = err;
        fprintf('Error 1: %f\n', err);

        % 2 
        fprintf('Running EM 2\n');
        [startTMine2, TMine2, EMine2, likelihood, gamma2] = EM2(X, m, n, 400);
        like2(end+1) = sum(log(gamma2(vec2mat(Y,2))))
        
        fprintf('Viterbi\n');
        YEst = viterbi(X, startTMine2, TMine2, EMine2);
        err = calcError(Y, YEst);
        e2(end + 1) = err;
        fprintf('Error 2: %f\n', err);


        % %1 WWW
        % fprintf('Running EM 3\n');
        % [model,log_like] = HMM_EM(X, m);
        
        % fprintf('Viterbi\n');
        % YEst = viterbi(X, model.P, model.A, model.B);
        % err = calcError(Y, YEst);
        % fprintf('Error 3: %f\n', err);
        % e3(end + 1) = err;

        plot(e1)
        hold on 
        plot(e2)
        % hold on 
        % plot(e3)
        hold off
        ylim([0,1]) ;
        legend('1st Order', '2nd Order');%, '1st Order from WWW');
        drawnow
    end
    figure
    E1(:,:) = EMine2(1,:,:);
    E2(:,:) = EMine2(2,:,:);
    subplot(1,2,1); imagesc(E1);
    title('Enhancer: est')
    subplot(1,2,2); imagesc(E2);
    title('Non-Enhancer: est')
    sum(EMine2(:))
    
    % figure
    % bar([mean(e1), mean(e2)]);%, mean(e3)])
    % title('EM 1st vs EM 2nd')

    figure;
    plot(like1);
    hold on;
    plot(like2);
    hold off
    legend('Log Likelihood 1st order', 'Log Likelihood 2nd order');
    title('Likelihood per repetition')
end


function XProb = getProb(X, n)
    k = length(X);
    XProb = sum(vec2mat(X, n), 2) ./ k;
    if k==0
        XProb = ones(n, 1) ./ n;
    end
end

function [startT, T, E] = genRandParams(m, n)
    [startT, T, E] = genRandParamsExt(m, n);
    [T, E] = addPrior(startT, T, E);
end

function [T1, E1] = addPrior(startT, T, E)
    [m, n] = size(E);
    E1 = [zeros(1, n); E];
    T1 = [0, startT.'; zeros(m, 1), T];
end


% E - 2 x 4 x 4
% T - 2 x 2
% startT- 1 x 2
function [startT, T, E] = genGenomeParams(epsilonT)
    % enhancer - 1
    % A C G T - 1 2 3 4
    m = 2;
    E1 = [223221, 156954, 239372, 156266; ...
          231040, 208677, 47147,  238196; ...
          190830, 169836, 208624, 156768; ...
          130811, 190164, 231375, 220719] / 3000000;
         % AA AC AG AT
         % CA CC CG CT
         % GA GC GG GT 
         % TA TC TG TT
    % non enhancer - 2
    E2= [300928, 158555, 210655, 252209; ...
         220123, 133040, 17280,  209043; ...
         182208, 107488, 131800, 157365; ...
         218356, 181179, 220291, 299480] / 3000000;
    E1 = bsxfun(@times, E1, 1 ./ sum(E1, 2));
    E2 = bsxfun(@times, E2, 1 ./ sum(E2, 2));
    E = permute(cat(3, E1, E2), [3,1,2]);

    % normalized random probabilities
    startT = rand(m, 1);
    startT = startT / sum(startT);
    % m x m 
    % T = eye(m);
    T = eye(m) + epsilonT * rand(m);
    T = bsxfun(@times, T, 1 ./ sum(T, 2));
    subplot(1,2,1); imagesc(E1);
    title('Enhancer: real')
    subplot(1,2,2); imagesc(E2);
    title('Non-Enhancer: real')
    sum(E(:))
    figure
end

% startT- 1 x 2
function [startT, T, E] = genStrong2dParams(epsilonT)
    % enhancer - 1
    % A C G T - 1 2 3 4
    m = 2;
    L = 30;
    E1 = [1, L, 1, 1; ...
          1, 1, L, 1; ...
          1, 1, 1, L; ...
          L, 1, 1, 1];
         % AA AC AG AT
         % CA CC CG CT
         % GA GC GG GT 
         % TA TC TG TT
    % non enhancer - 2
    E2 = [1, 1, 1, L; ...
          L, 1, 1, 1; ...
          1, L, 1, 1; ...
          1, 1, L, 1];
    
    E1 = bsxfun(@times, E1, 1 ./ sum(E1, 2));
    E2 = bsxfun(@times, E2, 1 ./ sum(E2, 2));
    E = permute(cat(3, E1, E2), [3,1,2]);

    % normalized random probabilities
    startT = rand(m, 1);
    startT = startT / sum(startT);
    % m x m 
    % T = eye(m);
    T = eye(m) + epsilonT * rand(m);
    T = bsxfun(@times, T, 1 ./ sum(T, 2));
    subplot(1,2,1); imagesc(E1);
    title('Enhancer: real')
    subplot(1,2,2); imagesc(E2);
    title('Non-Enhancer: real')
    sum(E(:))
    figure
end


function M3d = mat2mat3d(M, n)
    [r,c] = size(M);
    I = repmat(permute([1:n], [1,3,2]), [r, c, 1]);
    M3d = double(I == repmat(M, [1, 1, n]));
end

% function sdMat = rowsSD(E1, E2)
%     n = size(E1, 1);
%     m = size(E2, 1);
%     sdMat = zeros(n, m);
%     for i=1:n
%         for j=1:m
%             sdMat(i,j) = norm(E1(i, :) - E2(j, :));
%         end
%     end
%     figure
%     subplot(1,3,1);imagesc(sdMat);
%     title('Rows Match')
%     subplot(1,3,2);imagesc(E1);
%     title('Real E')
%     subplot(1,3,3);imagesc(E2);
%     title('Est. E')
%     drawnow
% end

% k - length of emmision 
% T - m x m transfer matrix
% E - m x n x n emission matrix
function [X, Y, likelihood] = genHmm2nd(k, startT, T, E)
    Y = zeros(1, k);
    X = zeros(1, k);
    Y(1) = find(mnrnd(1, startT));
    EStart(1, :) = sum(E(Y(1), :, :), 2);
    EStart = EStart / sum(EStart);
    X(1) = find(mnrnd(1, EStart));
    for i = 2:k
        Y(i) = find(mnrnd(1, T(Y(i-1), :)));
        ECur(:, 1) = E(Y(i), X(i-1), :);
        X(i) = find(mnrnd(1, ECur));
    end
end

