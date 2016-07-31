% make an EM version that uses context
% using real data of enhancer and non enhancer statistics


% Steps:
% get E T context matrix
% generate text
% EM guess the best E T 

% first version:
% get 1st level E T
% viterbi - get states Y
% get couples statistics, build 3D E T

% second version:
% really change EM file
function likelihoodPerM = ex4()
    close all;
    m = 2;
    n = 4;
    epsilonT = 0.05;

    k = 20000;
    fprintf('Samples length: %f\n', k);
    [startTReal, TReal, EReal] = [startT, T, E] = genGenomeParams(epsilonT);
    % Create X and Y
    [X, Y] = genHmm2nd(k, startTReal, TReal, EReal);
    [startTMine, TMine, EMine, likelihood] = EM(X, ms(i), n, 400);

    % estimate y
    YEst = viterbi(X, startTMine, TMine, EMine);

    % grade from -1 to 1 about how good is YEst
    % likelihood(end + 1) = silhouette(X, YEst, startTMine, TMine, EMine);
    dataPerM(i, t, 1) = likelihood;
    dataPerM(i, t, 2) = chi2(YEst, X, n, k, ms(i));
    dataPerM(i, t, 3) = gTest(YEst, X, n, k, ms(i));

        
    end

    likelihoodPerM(1,:) = mean(dataPerM(:, :, 1), 2)
    chi2PerM(1,:) = mean(dataPerM(:, :, 2), 2)
    gTestPerM(1,:) = mean(dataPerM(:, :, 3), 2)

    % log likelihood plot
    figure;
    plot(ms(1:i), likelihoodPerM(1:i),'b',[mReal], ...
        likelihoodPerM(mRealInd),'b*');
    title('Log(Likelihood) vs #states');
    xlabel('#states');
    ylabel('Log(Likelihood)');

    % model grades
    BITS_PER_NUMBER = 2^5;
    paramsAmounts = (ms .* (ms - 1) + ms .* (n - 1) + (ms - 1)) .* BITS_PER_NUMBER;
    bic = -2 .* likelihoodPerM + paramsAmounts .* log(k);
    aic = -2 .* likelihoodPerM + 2 .* paramsAmounts;
    figure;
    plot(ms, bic);
    hold on;
    plot(ms, aic);
    hold off;
    line([mReal mReal], get(gca,'YLim'), 'Color','r');
    title('Models Quality I');
    legend('BIC', 'AIC')
    
    figure;
    plot(ms, chi2PerM);
    hold on;
    plot(ms, gTestPerM);
    hold off;
    line([mReal mReal], get(gca,'YLim'), 'Color','r');
    title('Models Quality II');
    legend('chi2', 'gTest')

    LRT = -2.*(likelihoodPerM(1:end-1) - likelihoodPerM(2:end))
    figure;
    plot(ms(2:end), LRT);
    hold off;
    line([mReal mReal], get(gca,'YLim'), 'Color','r');
    title('Models Quality III');

    % save('C:\Users\booga\Dropbox\SemesterF\genetics\ex1\x.mat')
end
function [E, O] = chi2Ox(YEst, X, n, k, m)
    O = zeros(n, m);
    Xprob(:, 1) = getProb(X, n);
    E = repmat(Xprob, [1, m]);
    for i = 1 : m
        O(:, i) = getProb(X(YEst == i), n);
    end
end
function ret = gTest(YEst, X, n, k, m)
    [E, O] = chi2Ox(YEst, X, n, k, m);
    relInd = O~=0;
    ret = 2 * sum(O(relInd) .* log(O(relInd) ./ E(relInd)));
end

function ret = chi2(YEst, X, n, k, m)
    [E, O] = chi2Ox(YEst, X, n, k, m);
    relInd = E~=0;
    ret = ((O(relInd) - E(relInd)) .^ 2) ./ E(relInd);
    ret = sum(ret(:));
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
    % non enhancer - 2
    % A C G T - 1 2 3 4
    E1 = [223221, 156954, 239372, 156266; ...
         231040, 208677, 47147, 238196; ...
         190830, 169836, 208624, 156768; ...
         130811, 190164, 231375, 220719] / 3000000;
         % AA AC AG AT
         % CA CC CG CT
         % GA GC GG GT 
         % TA TC TG TT
    E2= [300928, 158555, 210655, 252209; ...
         220123, 133040, 17280, 209043; ...
         182208, 107488, 131800, 157365; ...
         218356, 181179, 220291, 299480] / 3000000;
    E = permute(cat(3, E1, E2), [3,1,2]);

    % normalized random probabilities
    startT = rand(m, 1);
    startT = startT / sum(startT);
    % m x m 
    % T = eye(m);
    T = eye(m) + epsilonT * rand(m);
    T = bsxfun(@times, T, 1 ./ sum(T, 2));
end


function M3d = mat2mat3d(M, n)
    [r,c] = size(M);
    I = repmat(permute([1:n], [1,3,2]), [r, c, 1]);
    M3d = double(I == repmat(M, [1, 1, n]));
end

function sdMat = rowsSD(E1, E2)
    n = size(E1, 1);
    m = size(E2, 1);
    sdMat = zeros(n, m);
    for i=1:n
        for j=1:m
            sdMat(i,j) = norm(E1(i, :) - E2(j, :));
        end
    end
    figure
    subplot(1,3,1);imagesc(sdMat);
    title('Rows Match')
    subplot(1,3,2);imagesc(E1);
    title('Real E')
    subplot(1,3,3);imagesc(E2);
    title('Est. E')
    drawnow
end

% k - length of emmision 
% T - m x m transfer matrix
% E - m x n x n emission matrix
function [X, Y, likelihood] = genHmm2nd(k, startT, T, E)
    Y = zero(1, k);
    X = zero(1, k);
    Y(1) = find(mnrnd(1, startT));
    EStart(1, :) = sum(E(Y(1), :, :), 2);
    EStart = EStart / sum(EStart)
    X(1) = find(mnrnd(1, EStart));
    for i = 2:k
        Y(i) = find(mnrnd(1, T(Y(i-1), :)));
        ECur(:, 1) = E(Y(i-1), X(i-1), :);
        ECur = ECur / sum(ECur);
        X(i) = find(mnrnd(1, ECur));
    end
end

