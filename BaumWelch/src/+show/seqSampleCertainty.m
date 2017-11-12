function seqSampleCertainty(params, Y, certainty)
    [N, L] = size(Y);
    sequencesToShow = 10;
    myColormap = repmat([1, 0.5, 0.5], [params.m, 1]);
    myColormap = [jet(params.m); myColormap];
    inds = randsample(N, sequencesToShow);
    inds = sort(inds);
    figure;
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
            title(['Posterior of ', num2str(sequencesToShow), ' Sequences']);
        end
        if i == sequencesToShow
            xlabel('Position in Sequence');
        else
            set(gca,'xtick',[]);
        end
    end
    legend(strcat({'Tissue Type '}, num2str([1:params.m]')))
end
function rotateYLabel()
    ylh = get(gca,'ylabel');
    gyl = get(ylh);
    ylp = get(ylh, 'Position');
    set(ylh, 'Rotation',0, 'Position',ylp, 'VerticalAlignment','middle', 'HorizontalAlignment','right')
end