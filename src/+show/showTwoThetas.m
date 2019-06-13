function showTwoThetas(params, thetaOrig, thetaEst, withExponent, subtitle, outpath)
    DOT_SIZE = 20;
    thetaEst = permThetaByAnother(params, thetaOrig, thetaEst);
    thetaOrigMat = thetaToMat(params, thetaOrig);
    thetaEstMat = thetaToMat(params, thetaEst);
    inds = [params.m, params.n ^ params.order, params.k, 1];
    colors = ['b', 'r', 'g', 'm'];
    fig = figure;
    if withExponent
        thetaOrigMat = exp(thetaOrigMat);
        thetaEstMat = exp(thetaEstMat);
    end
    minVal = floor(min([thetaEstMat(:); thetaOrigMat(:)]));
    maxVal = ceil(max([thetaEstMat(:); thetaOrigMat(:)]));
    plot([minVal, maxVal], [minVal, maxVal]);
    hold on;
    for i=1:4
        origMat = thetaOrigMat(:, 1:inds(i));
        thetaOrigMat = thetaOrigMat(:, inds(i) + 1:end);
        estMat = thetaEstMat(:, 1:inds(i));
        thetaEstMat = thetaEstMat(:,  inds(i) + 1:end);
        scatter(origMat(:), estMat(:), DOT_SIZE, colors(i), 'filled');
    end
    legend('x=y', 'T', 'E', 'G', 'startT', 'Location', 'southeast');
    title(sprintf('Learned \theta vs True \theta (%s)', subtitle));
    xlabel('True \theta');
    ylabel('Estimated \theta');
    saveas(fig, outpath);
end


function mat = thetaToMat(params, theta)
    mat = zeros(params.m, 1 + params.m + (params.n ^ params.order) + params.k);
    for i = 1:params.m
        mat(i, :) = [theta.T(i, :), theta.E(i, :), theta.G(i, :), theta.startT(i)];
    end
end


function theta = permThetaByAnother(params, thetaOrig, thetaEst)
    perm = findCorrectThetaPermute(params, thetaOrig, thetaEst);
    theta = permTheta(thetaEst, perm);
end

function theta = permTheta(theta, perm)
    theta.T = theta.T(perm, :);
    theta.T = theta.T(:, perm);
    theta.startT = theta.startT(perm);
    theta.G = theta.G(perm, :);
    for i = 1:length(perm)
        theta.E(i, :) = theta.E(perm(i), :);
    end
end

% perm - m x 1
function perm = findCorrectThetaPermute(params, thetaOrig, thetaEst)
    % to vec
    vectorizedOrig = thetaToMat(params, thetaOrig);
    vectorizedEst = thetaToMat(params, thetaEst);
    distMat = vectorizedOrig * vectorizedEst';
    % distMat = squareform(pdist(vectorized));
    perm = zeros(params.m, 1);
    for i = 1:params.m
        [~, I] = max(distMat(:));
        [I_row, I_col] = ind2sub(size(distMat), I);
        perm(I_col) = I_row;
        distMat(I_row, :) = -1;
        distMat(:, I_col) = -1;
    end
end

