% certainty - N x L
% Y - N x L
%posterior - N x m x L
function seqSampleCertainty(params, Y, gamma, psi)
    [N, L] = size(Y);
    % N x m x L
    posterior = calcPosterior(params, gamma, psi);
    % N x m x L
    YOneHot = permute(matUtils.mat23Dmat(Y, params.m), [1, 3, 2]);
    % N x L
    certainty = permute(sum(posterior .* YOneHot, 2), [1,3,2]);
    loss = mean(log(1-certainty(:)+eps));
    fprintf('Avg log loss: %.2f\n', loss);
    sequencesToShow = 10;
    myColormap = repmat([1, 0.5, 0.5], [params.m, 1]);
    myColormap = [jet(params.m); myColormap];
    inds = randsample(N, sequencesToShow);
    inds = sort(inds);
    loss = mean(log(1-certainty(:)));
    figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:sequencesToShow
        subplot(sequencesToShow, 1, i);
        YOneHot = matUtils.vec2mat(Y(inds(i), :), params.m);
        repCertainty = YOneHot .* repmat(certainty(inds(i), :), [params.m, 1]);
        repCertainty = [repCertainty; (1-repCertainty) .* double(repCertainty~=0)];
        bar(repCertainty', 'stacked')
        ylabel(['Seq ', num2str(inds(i))]);
        rotateYLabel();
        xlim([1, L])
        colormap(myColormap);
        if i == 1
            title(['Posterior of ', num2str(sequencesToShow), ' Sequences (logloss: ', num2str(loss),')']);
        end
        if i == sequencesToShow
            xlabel('Position in Sequence');
        else
            set(gca,'xtick',[]);
        end
    end
    legend(strcat({'Tissue Type '}, num2str([1:params.m]')));
    drawnow;
end

function rotateYLabel()
    ylh = get(gca,'ylabel');
    gyl = get(ylh);
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
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
