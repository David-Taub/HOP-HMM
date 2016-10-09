function addTissuesTicks(withAll, withBackground, ord, withY)
    tissues = {'All', 'BAT', 'BMDM', 'BoneMarrow',...
           'CH12', 'Cerebellum', 'Cortex',...
           'E14', 'Heart-E14.5', 'Heart',...
           'Kidney', 'Limb-E14.5', 'Liver-E14.5',...
           'Liver', 'MEF', 'MEL', 'OlfactBulb',...
           'Placenta', 'SmIntestine', 'Spleen',...
           'Testis', 'Thymus', 'WholeBrain-E14.5', 'mESC', 'background'};
    if ~withAll
        tissues = tissues(2:end);
    end
    if ~withBackground
        tissues = tissues(1:end-1);
    end
    M = length(tissues);
    if ord == 0
        ord = 1:M;
    end
    ax = gca;
    ax.XTick = 1:M;
    ax.XTickLabel = tissues(ord);
    ax.XTickLabelRotation=45;
    if withY
        ax.YTick = 1:M;
        ax.YTickLabel = tissues(ord);
    end
end
