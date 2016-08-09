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
    repeat = 3; 
    maxIter = 200; 
    tEpsilon = 0.005;
    m = 2;
    n = 4;
    k = 200000;
    
    e1 = zeros(1, repeat); 
    e2 = zeros(1, repeat);
    like1 = zeros(1, repeat); 
    like2 = zeros(1, repeat);
    
    [startTReal, TReal, EReal] = genGenomeParams(tEpsilon);
    % [startTReal, TReal, EReal] = genStrong2dParams(tEpsilon);
    % [startTReal, TReal, EReal] = genRandEyeParamsExt(m, n, tEpsilon, tEpsilon);

    for repIndex =1:repeat

        % Create X and Y
        fprintf('Genarated Params\n');
        % [startTReal, TReal, EReal]
        fprintf('Creating X and Y length: %d\n', k);
        [X, Y] = genHmm2nd(k, startTReal, TReal, EReal);
        % [X, Y] = genHmm(k, startTReal, TReal, EReal);

        % 1
        fprintf('Running EM\n');
        tic
        [startTEst1, TEst1, EEst1, likelihood, gamma] = EM(X, m, n, maxIter, tEpsilon);
        toc
        [like1(repIndex) , perm] = getLikelihood(Y, gamma);
        [startTEst1, TEst1, EEst1] = permuteParams(perm, startTEst1, TEst1, EEst1);

        fprintf('EM 1: %f\n', like1(repIndex));
        YEst1 = viterbi(X, startTEst1, TEst1, EEst1);
        err = calcError(Y, YEst1);
        e1(repIndex) = err;
        fprintf('Error 1: %f\n', err);

        % 2 
        fprintf('Running EM2\n');
        tic
        [startTEst2, TEst2, EEst2, likelihood, gamma2] = EM2(X, m, n, maxIter, tEpsilon);
        toc
        [like2(repIndex), perm] = getLikelihood(Y, gamma2);
        [startTEst2, TEst2, EEst2] = permuteParams(perm, startTEst2, TEst2, EEst2);
        fprintf('EM 2: %f\n', like2(repIndex));
        
        YEst2 = viterbi(X, startTEst2, TEst2, EEst2);
        err = calcError(Y, YEst2);
        e2(repIndex) = err;
        fprintf('Error 2: %f\n', err);
        TReal
        TEst1
        TEst2
        plotY(Y, YEst1, YEst2);
        plotE(EReal, EEst2);
    end
    plotLike(like1, like2, k);
    plotErr(e1, e2);
end
function plotY(Y, YEst1, YEst2)
    figure;
    hold on
    v = genDiffScatter(Y, 1);
    scatter(v, ones(1, length(v)) * 0.1, ones(1,length(v)) * 10, 's', 'filled');
    v = genDiffScatter(abs(diff(Y)), 0);
    scatter(repmat(v, [5, 1]), repmat([0.1: 0.1: 0.5].', [1, length(v)]), ones(1,length(v)) * 100, '+');
    v = genDiffScatter(YEst1, 1);
    scatter(v, ones(1, length(v)) * 0.2, ones(1,length(v)) * 10, 's', 'filled');
    v = genDiffScatter(YEst2, 1);
    scatter(v, ones(1, length(v)) * 0.3, ones(1,length(v)) * 10, 's', 'filled');
    v = genDiffScatter(YEst2, Y);
    scatter(v, ones(1, length(v)) * 0.4, ones(1,length(v)) * 10, 's', 'filled');
    v = genDiffScatter(YEst1, Y);
    scatter(v, ones(1, length(v)) * 0.5, ones(1,length(v)) * 10, 's', 'filled');
    hold off;
    ylim([0,1]);
    legend('Y', 'change Y', 'YEst1', 'YEst2', 'YEst1 Errors', 'YEst2 Errors');
end
function plotLike(like1, like2, k)
    figure
    plot(like1);
    hold on;
    plot(like2);
    legend('Likelihood GeoAvg 1', 'Likelihood GeoAvg 2');
    title('Likelihood GeoAvg per Repetition');
    hold off
end

function plotErr(e1, e2)
    figure
    plot(e1);
    hold on;
    plot(e2);
    legend('Vit Error 1', 'Vit Error 2');
    title('Viterbi Error per Repetition');
    hold off
end

function plotE(EReal, EEst2)
    figure
    E1(:,:) = EReal(1,:,:);
    E2(:,:) = EReal(2,:,:);
    E3(:,:) = EEst2(1,:,:);
    E4(:,:) = EEst2(2,:,:);
    subplot(2,2,1); imagesc(E1);
    title('Enhancer: real')
    subplot(2,2,2); imagesc(E2);
    title('Non-Enhancer: real')
    subplot(2,2,3); imagesc(E3);
    title('Enhancer: est')
    subplot(2,2,4); imagesc(E4);
    title('Non-Enhancer: est')
end

function v = genDiffScatter(Y1, Y2)
    L = 1; r = 1:length(Y1); v = r(Y1 ~= Y2) ./ L;
end

function [likelihood, perm]= getLikelihood(Y, gamma)
    like1 = sum(log(gamma(vec2mat(Y,2))));
    Y(Y==2) = 3;Y(Y==1) = 2; Y(Y==3) = 1;
    like2 = sum(log(gamma(vec2mat(Y,2))));
    if like1 > like2
        likelihood = like1;
        perm = [1, 2];
    else
        likelihood = like2;
        perm = [2, 1];
    end
    likelihood = exp(likelihood ./ length(Y));
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
function [startT, T, E] = genGenomeParams(tEpsilon)
    % enhancer - 1
    % A C G T - 1 2 3 4
    m = 2;
    E1 = [223221, 156954, 239372, 156266; ...
          231040, 208677, 47147,  238196; ...
          190830, 169836, 208624, 156768; ...
          130811, 190164, 231375, 220719];
         % AA AC AG AT
         % CA CC CG CT
         % GA GC GG GT 
         % TA TC TG TT
    % non enhancer - 2
    E2= [300928, 158555, 210655, 252209; ...
         220123, 133040, 17280,  209043; ...
         182208, 107488, 131800, 157365; ...
         218356, 181179, 220291, 299480];
    E1 = bsxfun(@times, E1, 1 ./ sum(E1, 2));
    E2 = bsxfun(@times, E2, 1 ./ sum(E2, 2));
    E = permute(cat(3, E1, E2), [3,1,2]);

    % normalized random probabilities
    startT = rand(m, 1);
    startT = startT / sum(startT);
    % m x m 
    % T = eye(m);
    noise1 = (tEpsilon * rand(1) ) * [-1, 1];
    noise2 = (tEpsilon * rand(1) ) * [1, -1];
    T = eye(2) + [noise1 ; noise2];
    % T = bsxfun(@times, T, 1 ./ sum(T, 2));
end

% startT- 1 x 2
function [startT, T, E] = genStrong2dParams(tEpsilon)
    % enhancer - 1
    % A C G T - 1 2 3 4
    m = 2;
    L = 5;
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
    T = eye(m) + tEpsilon * rand(m);
    T = bsxfun(@times, T, 1 ./ sum(T, 2));
    % subplot(1,2,1); imagesc(E1);
    % title('Enhancer: real')
    % subplot(1,2,2); imagesc(E2);
    % title('Non-Enhancer: real')
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



function [startT, T, E] = genRandEyeParamsExt(m, n, tEpsilon, epsilonE)
    % normalized random probabilities
    startT = rand(m, 1);
    startT = startT / sum(startT);
    % m x m 
    % T = eye(m);
    T = eye(m) + tEpsilon * rand(m);
    T = bsxfun(@times, T, 1 ./ sum(T, 2));
    % m x n
    E = eye(m, n) + epsilonE * rand(m, n);
    E = bsxfun(@times, E, 1 ./ sum(E, 2));
end




function [permStartT, permT, permE] = permuteParams(perm, startT, T, E)
    for i = 1 : length(perm)
        if length(size(E)) == 2
            permE(perm(i), :) = E(i, :);
        else
            permE(perm(i), :, :) = E(i, :, :);
        end
        permStartT(perm(i), 1) = startT(i);
        for j = 1 : length(perm)
            permT(perm(i), perm(j)) = T(i, j);
        end
    end
end