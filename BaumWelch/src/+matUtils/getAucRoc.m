
function auc = getAucRoc(pos, neg, shouldPlot)
    scores = [pos;neg];
    labels = [ones(length(pos), 1); zeros(length(neg), 1)];
    [X, Y, ~, auc] = perfcurve(labels, scores, 1);
    if shouldPlot
        plot(X,Y)
        title(['ROC AUC = ', sprintf('%.2f', auc)])
        xlabel('False positive rate')
        ylabel('True positive rate')
    end
end
