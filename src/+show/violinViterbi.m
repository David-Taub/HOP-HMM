function violinViterbi(truePosVals, trueNegVals, estPosVals, estNegVals, outpath)
    fig = figure();
    maxRange = max([truePosVals; trueNegVals; estPosVals; estNegVals], [], 1);
    minRange = min([truePosVals; trueNegVals; estPosVals; estNegVals], [], 1);
    show.distributionPlot({[truePosVals; minRange; maxRange], [trueNegVals; minRange; maxRange]}, 'histOri', 'right', 'color', 'r', 'widthDiv', [2 2], 'showMM', 0, 'histOpt', 0, 'divFactor', [minRange: maxRange]);
    show.distributionPlot({[estPosVals; minRange; maxRange], [estNegVals; minRange; maxRange]}, 'histOri', 'left', 'color', 'b', 'widthDiv', [2 1], 'showMM', 0, 'histOpt', 0, 'divFactor', [minRange: maxRange]);
    jf = java.text.DecimalFormat;
    tick1 = sprintf('TFBS Positions [n=%s, n=%s]',  char(jf.format(length(truePosVals))),  char(jf.format(length(estPosVals))));
    tick2 = sprintf('Non-TFBS Positions [n=%s, n=%s]',  char(jf.format(length(trueNegVals))),  char(jf.format(length(estNegVals))));
    set(gca,'xtick',[1:2],'xticklabel', {tick1, tick2});
    ylabel({'Log Posterior';'Probability (log\gamma)'});
    ylh = get(gca,'ylabel');
    ylp = get(ylh, 'Position');
    ext=get(ylh,'Extent');
    set(ylh, 'Rotation',0, 'Position',ylp-[ext(3)/2 0 0], 'VerticalAlignment','middle', 'HorizontalAlignment','center')
    legend( {'True TFBS positions', 'Viterbi estimated TFBS positions'})
    title('Posterior Probability of TFBS and non-TFBS Locations');
    saveas(fig, outpath);
end