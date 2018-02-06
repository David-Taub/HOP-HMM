% pos - N1 x 1
% neg - N2 x 1
function [auc, wasSwitched, thresh] = getAucRoc(pos, neg, shouldPlot, shouldSwitch)
    wasSwitched = false;
    if isempty(pos) || isempty(neg)
        fprintf('Warning empty input. pos size %d, neg size %d', length(pos), length(neg))
        auc = 0;
        wasSwitched = false;
        keyboard
        return;
    end
    scores = [pos;neg];
    labels = [ones(length(pos), 1); zeros(length(neg), 1)];
    % [X, Y, T, auc, OPTROCPT] = perfcurve(labels, scores, 1);
    [X, Y, T, auc, OPTROCPT] = matUtils.perfcurve2(labels, scores, 1);
    if auc < 0.5 && shouldSwitch
        [pos, neg] = deal(neg, pos);
        scores = [pos;neg];
        labels = [ones(length(pos), 1); zeros(length(neg), 1)];
        wasSwitched = true;
        [X, Y, T, auc, OPTROCPT] = matUtils.perfcurve2(labels, scores, 1);
    end
    thresh = T((X==OPTROCPT(1))&(Y==OPTROCPT(2)));
    if shouldPlot
        figure
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
