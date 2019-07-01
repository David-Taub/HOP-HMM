% sample sequences, and draw for each colorful plots with what the posterior
% probability was compared to the correct state per letter
    function seqSampleCertaintyReal(params, theta, dataset, sequencesToShow, outpath)
    [N, L] = size(dataset.X);
    % gamma - N x m x L
    % psi - N x m x k x L
    [~, ~, ~, ~, gamma, psi] = EM.EStep(params, theta, dataset.X, dataset.pcPWMp);
    % N x m x L
    posterior = calcPosterior(params, gamma, psi);
    % sequencesToShow = 10;
    cMap = lines(params.m);
    cellCMap = {};
    for i=1:params.m
        cellCMap{i, 1} = cMap(i, :);
    end

    PWM_COLOR = [0, 0, 0];
    cMapWithError = [cMap;  PWM_COLOR];
    inds = randsample(N, sequencesToShow);
    inds = sort(inds);
    figure('units', 'pixels', 'Position', [0 0 1000 1000]);
    YEst = misc.viterbi(params, theta, dataset.X, dataset.pcPWMp);
    pwm_val = params.m + 1;
    LOW_BAR_HIEGHT = 0.1;
    for i = 1:sequencesToShow
        ind = inds(i);
        subplot(sequencesToShow, 1, i);
        hold on;
        % lowbar
        ylim([-LOW_BAR_HIEGHT, 1]);

        YEstOneHot = matUtils.vec2mat(YEst(ind, :, 1), params.m);
        YEst(ind, YEst(ind, :, 2) > 0, 1) = pwm_val;

        imagesc([.5, L-.5], [-3 * LOW_BAR_HIEGHT / 4, -LOW_BAR_HIEGHT / 4], [YEst(ind,:,1);YEst(ind,:,1)], [1, pwm_val]);
        colormap(cMapWithError);
        text(L + 1, 0.5, 'Posterior Probability', 'FontSize', 8)
        text(L + 1, -LOW_BAR_HIEGHT / 2, 'Viterbi Estimation', 'FontSize', 8)

        % m x L
        selectedPosterior = permute(posterior(ind, :, :), [2, 3, 1]) .* YEstOneHot;
        % m + 1 x L
        b = bar(1:L, selectedPosterior', 1, 'stacked', 'FaceColor','flat');
        ylabel(['Seq ', num2str(ind)]);

        rotateYLabel();
        xlim([1, L])
        for j = 1:size(selectedPosterior, 1)
            b(j).CData = cMap(j, :) * 0.85;
        end

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
