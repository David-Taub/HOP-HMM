% sample sequences, and draw for each colorful plots with what the posterior
% probability was compared to the correct state per letter
% Y - N x L
function seqSampleCertainty(params, theta, dataset, sequencesToShow, outpath)
    [N, L] = size(dataset.Y);
    % gamma - N x m x L
    % psi - N x m x k x L
    [~, ~, ~, ~, gamma, psi] = EM.EStep(params, theta, dataset.X, dataset.pcPWMp);
    % TODO: think about summing on the correct PWM only, not all of them
    % N x m x L
    posterior = calcPosterior(params, gamma, psi);
    % N x m x L
    YOneHot = permute(matUtils.mat23Dmat(dataset.Y, params.m), [1, 3, 2]);
    % N x L - the posteriors of the correct state
    certainty = permute(sum(posterior .* YOneHot, 2), [1,3,2]);
    loss = mean(1-certainty(:) + eps);
    fprintf('Avg loss: %.2f\n', loss);
    % sequencesToShow = 10;
    cMap = lines(params.m);
    cellCMap = {};
    for i=1:params.m
        cellCMap{i, 1} = cMap(i, :);
    end
    YEst = misc.viterbi(params, theta, dataset.X, dataset.pcPWMp);
    ERROR_COLOR = [1.0, 0.1, 0.1];
    PWM_COLOR = [0, 0, 0];
    MATCH_COLOR = [1, 1, 1];
    cMapWithError = [cMap; ERROR_COLOR; PWM_COLOR; MATCH_COLOR];
    inds = randsample(N, sequencesToShow);
    inds = sort(inds);
    loss = mean(1-certainty(:));
    figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:sequencesToShow
        ind = inds(i);
        subplot(sequencesToShow, 1, i);
        hold on;
        % lowbar
        match_val = params.m + 3;
        pwm_val = params.m + 2;
        error_val = params.m + 1;
        LOW_BAR_HIEGHT = 0.2;

        ylim([-LOW_BAR_HIEGHT, 1]);
        errorMask = any(YEst(ind, :, :) ~= cat(3, dataset.Y(ind, :), dataset.Y2(ind, :)), 3);
        errorHeat = errorMask .* error_val + not(errorMask) .* match_val;
        origY = dataset.Y(ind, :);
        origY(dataset.Y2(ind, :) > 0) = pwm_val;
        YEst(ind, YEst(ind, :, 2) > 0, 1) = pwm_val;
        heatmap = [origY; YEst(ind, :, 1)];
        imagesc([.5, L-.5], [-3 * LOW_BAR_HIEGHT / 4, -LOW_BAR_HIEGHT / 4], heatmap, [1, params.m + 3]);
        colormap(cMapWithError);

        im = imagesc([.5, L-.5], [-52 * LOW_BAR_HIEGHT / 100, -48 * LOW_BAR_HIEGHT / 100], [errorHeat; errorHeat], [1, params.m + 3]);
        colormap(cMapWithError);
        im.AlphaData = [errorMask;errorMask] > 0;
        text(L + 1, 0.5, 'Posterior Probability', 'FontSize', 8)
        text(L + 1, -LOW_BAR_HIEGHT / 4, 'Viterbi Estimation', 'FontSize', 8)
        text(L + 1, -3 * LOW_BAR_HIEGHT / 4, 'Real States', 'FontSize', 8)

        % heatmap = [origY; YEst(ind, :, 1); errorMask];
        % imagesc([.5, L-.5], [-5 * LOW_BAR_HIEGHT / 6, -LOW_BAR_HIEGHT / 6], heatmap, [1, params.m + 3]);

        % m x L
        YOneHot = matUtils.vec2mat(dataset.Y(ind, :), params.m);
        % m x L
        correctPosterior = permute(posterior(ind, :, :), [2, 3, 1]) .* YOneHot;
        incorrectPosterior = sum(permute(1 - posterior(ind, :, :), [2, 3, 1]) .* YOneHot, 1);
        % m + 1 x L
        barsInfo = [correctPosterior; incorrectPosterior];
        b = bar(1:L, barsInfo', 1, 'stacked', 'FaceColor','flat');
        ylabel(['Seq ', num2str(ind)]);

        rotateYLabel();
        xlim([1, L])
        for j = 1:size(barsInfo, 1)
            b(j).CData = cMapWithError(j, :) * 0.85;
        end
        b(params.m + 1).CData = ERROR_COLOR;

        if i == 1
            title('Posterior & Viterbi Estimation');
        end
        if i == sequencesToShow
            xlabel('Position in Sequence');
        else
            set(gca,'xtick',[]);
        end
        p = plot(1:L, permute(posterior(ind, :, :), [3, 2, 1]), 'LineWidth', 1.5);
        set(p, {'Color'}, cellCMap);
        yticks([0:0.2:1]);
        hold off;
    end
    legendStrings1 = strcat({'Enhancer Type '}, num2str([1:params.m - params.backgroundAmount]'));
    legendStrings2 = strcat({'Background '}, num2str([1:params.backgroundAmount]'));
    legendStrings = {legendStrings1{:}, legendStrings2{:}};
    legendStrings{error_val} = 'Error';
    legendStrings{pwm_val} = 'TFBS';
    legend(legendStrings);
    saveas(gcf, outpath);
end

function rotateYLabel()
    ylh = get(gca,'ylabel');
    gyl = get(ylh);
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation', 0, 'Position', ylp, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right')
end

% posterior - N x m x L
function posterior = calcPosterior(params, gamma, psi)
    N = size(gamma, 1);
    posterior = gamma;
    for l = 1:params.k
        for t = 1:params.lengths(l)
            subStatePost = permute(cat(4, -inf(N, params.m, 1, t), psi(:, :, l, 1:end-t)), [1,2,4,3]);
            posterior = matUtils.logAdd(posterior, subStatePost);
        end
    end
    posterior = exp(posterior);
end
