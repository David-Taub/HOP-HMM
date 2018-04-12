% certainty - N x L
% Y - N x L
%posterior - N x m x L
function seqSampleCertainty(params, Y, gamma, psi, sequencesToShow, showFirst)
    [N, L] = size(Y);
    % N x m x L
    posterior = calcPosterior(params, gamma, psi);
    % N x m x L
    YOneHot = permute(matUtils.mat23Dmat(Y, params.m), [1, 3, 2]);
    % N x L
    certainty = permute(sum(posterior .* YOneHot, 2), [1,3,2]);
    loss = mean(log(1-certainty(:)+eps));
    fprintf('Avg log loss: %.2f\n', loss);
    % sequencesToShow = 10;
    mJet = jet(params.m);
    cellJet = {};
    for i=1:params.m
        cellJet{i, 1} = mJet(i, :);
    end
    myColormap = [mJet; repmat([1, 0.5, 0.5], [params.m, 1])];
    inds = randsample(N, sequencesToShow);
    if showFirst
        inds = 1:sequencesToShow;
    end

    inds = sort(inds);
    loss = mean(log(1-certainty(:)));
    figure('units','normalized','outerposition',[0 0 1 1]);
    for i = 1:sequencesToShow
        subplot(sequencesToShow, 1, i);

        % m x L
        YOneHot = matUtils.vec2mat(Y(inds(i), :), params.m);
        % m x L
        repCertainty = YOneHot .* repmat(certainty(inds(i), :), [params.m, 1]);
        % 2*m x L
        repCertainty = [repCertainty; (1-repCertainty) .* double(repCertainty~=0)];
        bar(1:L, repCertainty', 1, 'stacked')
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
        hold on;
        p = plot(1:L, permute(posterior(inds(i), :, :), [3, 2, 1]));
        set(p, {'Color'}, cellJet)
        hold off;
    end
    legendStrings = strcat({'Tissue Type '}, num2str([1:params.m-1]'));
    legendStrings{params.m} = 'Background';
    legendStrings{params.m+1} = 'Error';
    legend(legendStrings);
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
