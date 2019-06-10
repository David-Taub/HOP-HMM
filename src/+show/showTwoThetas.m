function showTwoThetas(params, thetaOrig, thetaEst, withExponent)
    DOT_SIZE = 20
    thetaOrigMat = thetaToMat(params, thetaOrig);
    thetaEstMat = thetaToMat(params, thetaEst);
    inds = [params.m, params.n ^ params.order, params.k, 1]
    colors = ['b', 'r', 'g', 'm']
    figure;
    hold on;
    if withExponent
        thetaOrigMat = exp(thetaOrigMat);
        thetaEstMat = exp(thetaEstMat);
        minVal = floor(min([thetaEstMat(:); thetaOrigMat(:)]))
        maxVal = ceil(max([thetaEstMat(:); thetaOrigMat(:)]))
        plot([minVal, maxVal], [minVal, maxVal])
    end
    for i=1:4
        origMat = thetaOrigMat(:, 1:inds(i));
        thetaOrigMat = thetaOrigMat(:, inds(i)+1:end);
        estMat = thetaEstMat(:, 1:inds(i));
        thetaEstMat = thetaEstMat(:,  inds(i)+1:end);
        scatter(exp(origMat(:)), exp(estMat(:)), DOT_SIZE, colors(i), 'filled');
    end
    legend('x=y', 'T', 'E', 'G', 'startT')
    title('Learned Parameters vs True Parameters')
    xlabel('True')
    ylabel('Estimated')
end
