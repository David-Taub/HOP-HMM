% sample sequences, and draw for each colorful plots with what the posterior
% probability was compared to the correct state per letter
function seqSampleCertainty(params, theta, dataset, sequencesToShow, outpath)
    ERROR_COLOR = [1.0, 0.05, 0.05];
    PWM_COLOR = [0, 0, 0];
    MATCH_COLOR = [1, 1, 1];

    [N, L] = size(dataset.Y);
    % gamma - N x m x L
    % psi - N x m x k x L
    [~, ~, ~, ~, gamma, psi] = EM.EStep(params, theta, dataset.X, dataset.pcPWMp);
    % N x m x L
    TFPosterior = calcTFPosterior(params, psi);
    TFPosterior = matUtils.logMatSum(TFPosterior, 2);
    % N x m x L
    posterior = exp(cat(2, gamma, TFPosterior));
    % sequencesToShow = 10;
    colorMap = lines(params.m);
    cellColorMap = {};
    for i = 1:params.m
        cellColorMap{i, 1} = colorMap(i, :);
    end
    cellColorMap{params.m + 1, 1} = PWM_COLOR;
    % YEst - N x L x 2
    YEst = misc.viterbi(params, theta, dataset.X, dataset.pcPWMp);



    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    confusionOutPath = "confMat_" + outpath;
    Y = dataset.Y;
    Y(:, :, 2) = dataset.Y2;
    showConfusionMatrix(params, YEst, Y);
    saveas(gcf, confusionOutPath);



    colorMapWithError = [colorMap; PWM_COLOR; ERROR_COLOR; MATCH_COLOR];
    inds = randsample(N, sequencesToShow);
    inds = sort(inds);
    % lowbar
    pwmVal = params.m + 1;
    errorVal = params.m + 2;
    matchVal = params.m + 3;
    LOW_BAR_HIEGHT = 0.2;
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    for i = 1:sequencesToShow
        ind = inds(i);
        subplot(sequencesToShow, 1, i);
        hold on;

        % ############################
        % Plot 2 bars at bottom
        % ############################

        ylim([-LOW_BAR_HIEGHT, 1]);
        errorMask = any(YEst(ind, :, :) ~= cat(3, dataset.Y(ind, :), dataset.Y2(ind, :)), 3);
        errorHeat = errorMask .* errorVal + not(errorMask) .* matchVal;
        origY = dataset.Y(ind, :);
        origY(dataset.Y2(ind, :) > 0) = pwmVal;
        YEst(ind, YEst(ind, :, 2) > 0, 1) = pwmVal;
        heatmap = [origY; YEst(ind, :, 1)];
        imagesc([.5, L - .5], [-3 * LOW_BAR_HIEGHT / 4, -LOW_BAR_HIEGHT / 4], heatmap, [1, params.m + 3]);
        colormap(colorMapWithError);
        % separator line
        plot([0, L], [-LOW_BAR_HIEGHT / 2, -LOW_BAR_HIEGHT / 2], 'LineWidth', 1, 'Color', 'k');

        text(L + 5, 0.5, 'Posterior Probability', 'FontSize', 10);
        text(L + 5, -LOW_BAR_HIEGHT / 4, 'Viterbi Estimation', 'FontSize', 10);
        text(L + 5, -3 * LOW_BAR_HIEGHT / 4, 'Real States', 'FontSize', 10);


        % ############################
        % Plot upper posterior section
        % ############################

        % m x L
        YOneHot = matUtils.vec2mat(dataset.Y(ind, :), params.m);
        YOneHot(:, dataset.Y2(ind, :) > 0) = 0;
        % m + 1 x L
        YOneHot = cat(1, YOneHot, dataset.Y2(ind, :) > 0);

        % m + 1 x L
        correctPosterior = permute(posterior(ind, :, :), [2, 3, 1]) .* YOneHot;
        % m + 1 x L
        incorrectPosterior = sum(permute(1 - posterior(ind, :, :), [2, 3, 1]) .* YOneHot, 1);
        % m + 1 x L
        barsValues = [correctPosterior; incorrectPosterior];

        barPlotHandler = bar(1:L, barsValues', 1, 'stacked', 'FaceColor','flat');
        ylabel(['Seq ', num2str(ind)]);

        rotateYLabel();
        xlim([1, L])
        for j = 1:size(barsValues, 1)
            barPlotHandler(j).CData = colorMapWithError(j, :) * 0.85;
        end
        barPlotHandler(params.m + 2).CData = ERROR_COLOR;

        if i == 1
            title('Posterior & Viterbi Estimation');
        end
        if i == sequencesToShow
            xlabel('Position in Sequence');
        else
            set(gca,'xtick', []);
        end
        p = plot(1: L, permute(posterior(ind, :, :), [3, 2, 1]), 'LineWidth', 1.5);
        set(p, {'Color'}, cellColorMap);
        yticks([0: 0.2: 1]);
        hold off;
    end
    legendStrings1 = strcat({'Enhancer Type '}, num2str([1:params.enhancerAmount]'));
    legendStrings2 = strcat({'Background '}, num2str([1:params.backgroundAmount]'));
    legendStrings = {legendStrings1{:}, legendStrings2{:}};
    legendStrings{errorVal} = 'Error';
    legendStrings{pwmVal} = 'TFBS';
    legend(legendStrings);
    saveas(gcf, outpath);
end


function rotateYLabel()
    ylh = get(gca,'ylabel');
    gyl = get(ylh);
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation', 0, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
end

% gamma - N x m x L
% psi - N x m x k x L
% enhacnerPosterior - N x m x L
% TFPosterior - N x m x L
function TFPosterior = calcTFPosterior(params, psi)
    N = size(psi, 1);
    L = size(psi, 4);
    TFPosterior = -inf(N, params.m, L);
    for l = 1:params.k
        for t = 1:params.lengths(l)
            % N x m x L
            TFStatePost = permute(cat(4, -inf(N, params.m, 1, t), psi(:, :, l, 1:end - t)), [1, 2, 4, 3]);
            TFPosterior = matUtils.logAdd(TFPosterior, TFStatePost);
        end
    end
end

% YEst - N x L x 2
% Y - N x L x 2
function showConfusionMatrix(params, YEst, Y)
    % plot a confusion matrix
    enhancerStateY = Y(:, :, 1);
    enhancerStateYEst = YEst(:, :, 1);
    confusionRaw = confusionmat(enhancerStateY(:), enhancerStateYEst(:));
    confusion = confusionRaw ./ repmat(sum(confusionRaw, 2), [1, size(confusionRaw, 2)]);
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    imagesc(confusion); colorbar; title('Viterbi Confusion Matrix (floors only)');
    caxis([0 1]);
    ax = gca;
    ax.XTick = [1 : size(confusion, 1)];
    ax.YTick = [1 : size(confusion, 1)];
    xlabel('Estimated Enhancer States');
    ylabel('True Enhancer States');


    enhancerStateY = 100 * Y(:, :, 1) + Y(:, :, 2);
    enhancerStateYEst = 100 * YEst(:, :, 1) + YEst(:, :, 2);
    confusionRaw = confusionmat(enhancerStateY(:), enhancerStateYEst(:));
    confusion = confusionRaw ./ repmat(sum(confusionRaw, 2), [1, size(confusionRaw, 2)]);
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    imagesc(confusion); colorbar; title('Viterbi Confusion Matrix (including TF states)');
    caxis([0 1]);
    ax = gca;
    ax.XTick = [1 : size(confusion, 1)];
    ax.YTick = [1 : size(confusion, 1)];
    xlabel('Estimated Enhancer States');
    ylabel('True Enhancer States');

    confusion = logRaw(confusionRaw);
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    imagesc(confusion); colorbar; title('Viterbi Confusion Matrix (log)');
    ax = gca;
    ax.XTick = [1 : size(confusion, 1)];
    ax.YTick = [1 : size(confusion, 1)];
    xlabel('Estimated Enhancer States');
    ylabel('True Enhancer States');
end