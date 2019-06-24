function showTwoThetas(params, thetaOrig, thetaEst, withExponent, subtitle, outpath)
    DOT_SIZE = 20;
    thetaEst = misc.permThetaByAnother(params, thetaOrig, thetaEst);
    thetaOrigMat = misc.thetaToMat(params, thetaOrig, true);
    thetaEstMat = misc.thetaToMat(params, thetaEst, true);
    inds = [params.m, params.n ^ params.order, params.k, 1];
    colors = ['b', 'r', 'g', 'm'];
    if isempty(strfind(outpath, 'tmp'))
        fig = figure('units','normalized','outerposition',[0 0 1 1]);
    end
    if withExponent
        thetaOrigMat = exp(thetaOrigMat);
         thetaEstMat = exp(thetaEstMat);
    end
    minVal = floor(min([thetaEstMat(:); thetaOrigMat(:)]));
    maxVal = ceil(max([thetaEstMat(:); thetaOrigMat(:)]));
    plot([minVal, maxVal], [minVal, maxVal]);
    hold on;
    % mse = mean((thetaOrigMat(:) - thetaEstMat(:)) .^ 2, 1);
    msle = mean((log(exp(thetaOrigMat(:)) + 1) - log(exp(thetaEstMat(:)) + 1)) .^ 2, 1);
    for i=1:4
        % plot by color
        origSubMat = thetaOrigMat(:, 1:inds(i));
        estSubMat = thetaEstMat(:, 1:inds(i));

        scatter(origSubMat(:), estSubMat(:), DOT_SIZE, colors(i), 'filled');

        thetaOrigMat = thetaOrigMat(:, inds(i) + 1:end);
        thetaEstMat = thetaEstMat(:,  inds(i) + 1:end);
    end
    legend('x=y', 'T', 'E', 'G', 'startT', 'Location', 'southeast');
    title(sprintf('Learned \\theta vs True \\theta (%s)', subtitle));
    xlabel('True \theta');
    ylabel('Estimated \theta');
    text(minVal + 0.1 * (maxVal - minVal), minVal + 0.8 * (maxVal - minVal), sprintf('MSLE=%.2e', msle), 'FontSize', 14)
    saveas(gcf, outpath);
end


function showTwoThetasOverTime(params, thetaOrig, thetaEsts, withExponent, subtitle, outpath)
    DOT_SIZE = 20;
    % colors = ['b', 'r', 'g', 'm'];
    hold on;
    BASE_INTENSITY = 0.7;
    for j = 1:length(thetaEsts)
        thetaEst = misc.permThetaByAnother(params, thetaOrig, thetaEsts(j));
        % m x length(theta(:)) / m
        thetaOrigMat = misc.thetaToMat(params, thetaOrig, true);
        thetaEstMat = misc.thetaToMat(params, thetaEst, true);
        if isempty(strfind(outpath, 'tmp'))
            fig = figure('units','normalized','outerposition',[0 0 1 1]);
        end
        if withExponent
            thetaOrigMat = exp(thetaOrigMat);
            thetaEstMat = exp(thetaEstMat);
        end
        minVal = floor(min(minVal, min([thetaEstMat(:); thetaOrigMat(:)])));
        maxVal = ceil(max(maxVal, max([thetaEstMat(:); thetaOrigMat(:)])));
        % mse = mean((thetaOrigMat(:) - thetaEstMat(:)) .^ 2, 1);
        msle = mean((log(exp(thetaOrigMat(:)) + 1) - log(exp(thetaEstMat(:)) + 1)) .^ 2, 1);
        gradColor = BASE_INTENSITY * (1 - (j / length(thetaEsts)));
        color = [gradColor, gradColor, 1];
        scatter(origSubMat(:), estSubMat(:), DOT_SIZE, color, 'filled');
        if j > 1
            for t = 1: length(origSubMat(:))
            plot(origSubMat(:), estSubMat(:), DOT_SIZE, color, 'filled');
        end
        prevTheta = thetaEstMat

        end
    plot([minVal, maxVal], [minVal, maxVal]);
    end
    % legend('x=y', 'T', 'E', 'G', 'startT', 'Location', 'southeast');
    title(sprintf('Learned \\theta vs True \\theta (%s)', subtitle));
    xlabel('True \theta');
    ylabel('Estimated \theta');
    text(minVal + 0.1 * (maxVal - minVal), minVal + 0.8 * (maxVal - minVal), sprintf('MSLE=%.2e', msle), 'FontSize', 14)
    saveas(gcf, outpath);
end