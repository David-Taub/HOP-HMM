
function auc = getAucRoc(pos, neg, shouldPlot)
    if length(pos) == 0 || length(neg) == 0
        auc = 0;
        return
    end

    if mean(pos, 1) < mean(neg, 1)
        [pos, neg] = deal(neg, pos);
    end
    scores = [pos;neg];
    labels = [ones(length(pos), 1); zeros(length(neg), 1)];
    [X, Y, ~, auc] = perfcurve(labels, scores, 1);
    if shouldPlot
        subplot(1,2,2);
        plot(X,Y)
        title(['ROC AUC = ', sprintf('%.2f', auc)])
        xlabel('False positive rate')
        ylabel('True positive rate')

        subplot(1,2,1);
        h = histogram(pos, 50, 'Normalization', 'probability');
        hold on;
        histogram(neg, h.BinEdges, 'Normalization', 'probability');
        legend('Enhancers of cell 1', 'Enhancers of cell 2')
        hold off;
        drawnow
    end
end
