
function selectedPWMs = PWMsFeatureSelect(mergedPeaksMin, wantedK)
    fprintf('Selecting PWMs\n');


    dbstop if error
    % PWM - k x n x J
    % lengths - k x 1
    % names - k x 1
    [PWMs, lengths, names] = misc.PWMs();
    m = size(mergedPeaksMin.overlaps, 2);
    L = size(mergedPeaksMin.seqs, 2);
    [k, n, J] = size(PWMs);
    outFilepath = sprintf('../data/precomputation/SelectedPWMs_L%dwk%dm%dtissues%s.mat', L,...
                          wantedK, m, sprintf('%dvs%d', mergedPeaksMin.tissueList(1), mergedPeaksMin.tissueList(2)));

    if isfile(outFilepath)
        fprintf('Found selected PWMS cache file at %s\n', outFilepath);
        load(outFilepath);
        return;
    end
    [X, Y] = preprocess(mergedPeaksMin);
    Xs1H = matUtils.mat23Dmat(X, n);


    selectedPWMs = [];
    selectedPWMs = getExpressedPWMs(mergedPeaksMin, PWMs, names, selectedPWMs, wantedK);
    fprintf('Selected PWMs by expression: %s\n', sprintf('%d ', selectedPWMs));
    [names{selectedPWMs}]
    if length(selectedPWMs) < wantedK
        selectedPWMs = getStrongPWMs(Xs1H, PWMs, names, lengths, wantedK, selectedPWMs);
        fprintf('Selected PWMs by strength: %s\n', sprintf('%d ', selectedPWMs));
        [names{selectedPWMs}]
        if length(selectedPWMs) < wantedK
            selectedPWMs = getDistinguishablePWMs(Xs1H, Y, PWMs, lengths, wantedK, selectedPWMs);
            fprintf('Added PWMs by AUC-ROC between sequences: %s\n', sprintf('%d ', selectedPWMs));
            [names{selectedPWMs}]
        end
    end
    save(outFilepath, 'selectedPWMs');
    fprintf('Saved feature selected PWMs in %s\n', outFilepath);
end

function selectedPWMs = getExpressedPWMs(mergedPeaksMin, PWMs, names, selectedPWMs, wantedK)
    PICK_EXPRESSED_RATIO = 0.3;
    MAX_REVERSE_MAP = 4;
    expressionTable = readtable('..\data\Jaspar\raw\TF_expressions.txt');
    TFNames = expressionTable.gene_id;
    PWMsNames = upper([names{:}])';

    % match pwms and TF expression names
    reverseIds = {};
    for i = 1:length(TFNames)
        ids = find(cellfun(@isempty, strfind(PWMsNames, TFNames{i})) == 0);
        if ~isempty(ids)
            if length(ids) > MAX_REVERSE_MAP
                ids = ids(1:MAX_REVERSE_MAP);
            end
            reverseIds{i} = ids';
        end
    end
    % find over expressed TFs
    eid1 = mergedPeaksMin.tissueEIDs{1};
    eid2 = mergedPeaksMin.tissueEIDs{2};
    eidMask1 = strcmp(eid1, expressionTable.Properties.VariableNames);
    eidMask2 = strcmp(eid2, expressionTable.Properties.VariableNames);
    if not(sum(eidMask1) == 1 & sum(eidMask2) == 1)
        fprintf('No expressed data for selected tissues\n');
        return
    end
    expression1 = expressionTable{:, eidMask1};
    psaudoCount = min(expression1(expression1 > 0));
    expression1 = (expression1 + psaudoCount) ./ norm(expression1);
    expression2 = expressionTable{:, eidMask2};
    expression2 = (expression2 + psaudoCount) ./ norm(expression2);
    tfExpressionRatio = expression1 ./ expression2;

    [~, inds] = sort(tfExpressionRatio, 'descend');
    expressedPWMsIds1 = [reverseIds{inds(1: floor(wantedK * PICK_EXPRESSED_RATIO / 2))}];
    expressedPWMsIds2 = [reverseIds{inds(end - floor(wantedK * PICK_EXPRESSED_RATIO / 2): end)}];
    selectedPWMs = [selectedPWMs, expressedPWMsIds1, expressedPWMsIds2];
end

% Xs1H - N x L x n
% Y - N x L
% PWMs - k x n x J
% lengths - 1 x J
function selectedPWMs = getDistinguishablePWMs2(Xs1H, Y, PWMs, lengths, wantedK, selectedPWMs)
    % aucRocs  - m x k
    aucRocs = misc.oneVsAllAucRoc(Xs1H, Y, PWMs, lengths);
    [aucRocsSorted, aucRocsSortedInd] = sort(aucRocs, 2, 'descend');
    i = 1;
    while length(selectedPWMs) < wantedK
        selectedPWMs = [selectedPWMs, unique(aucRocsSortedInd(1:i))];
        i = i + 1;
    end
end

% Xs1H - N x L x n
% Y - N x 1
% PWMs - k x n x J
% lengths - 1 x J
function selectedPWMs = getDistinguishablePWMs(Xs1H, Y, PWMs, lengths, wantedK, selectedPWMs)
    % aucRocs  - m x k
    [N, L, ~] = size(Xs1H);
    Xs1HNew = Xs1H(:, floor(L / 3) + 1: 2 * floor(L / 3), :);
    Xs1HNew = cat(1, Xs1HNew, Xs1H(:, 1: floor(L / 3), :));
    Xs1HNew = cat(1, Xs1HNew, Xs1H(:, 2 * floor(L / 3) + 1: 3 * floor(L / 3), :));
    m = max(Y(:));
    Y = cat(1, Y, (m + 1) .* ones(2 * N, 1));
    aucRocs = misc.oneVsAllAucRoc(Xs1HNew, Y, PWMs, lengths);
    [aucRocsSorted, aucRocsSortedInd] = sort(aucRocs, 2, 'descend');
    i = 1;
    while length(selectedPWMs) < wantedK
        selectedPWMs = [selectedPWMs, unique(aucRocsSortedInd(1:i))];
        i = i + 1;
    end
end

function selectedPWMs = getStrongPWMs(Xs1H, PWMs, names, lengths, wantedK, selectedPWMs)
    EXPECTED_NUM_OF_PEAKS_IN_SEQ = 3;
    PICK_STRONGEST_RATIO = 0.3;
    [k, n, ~] = size(PWMs);
    rates = zeros(1, k);
    for pwmId = 1:k
        PWMLogLike = misc.PWMLogLikelihood(PWMs, lengths, Xs1H, pwmId);
        % mask x EXPECTED_NUM_OF_PEAKS_IN_SEQ
        rates(pwmId) = mean(mean(misc.maxN(PWMLogLike, 2, EXPECTED_NUM_OF_PEAKS_IN_SEQ), 2), 1);
        rates(pwmId) = rates(pwmId) - log(1 / n) * lengths(pwmId);
        name = names{pwmId};
        fprintf('Checked PWM %d %s- log likelihood: %.2f (%.2f), length: %d\n', pwmId, name{:}, rates(pwmId), max(rates), lengths(pwmId));
    end
    [rates, inds] = sort(rates, 'descend');
    selectedPWMs = [selectedPWMs, inds(1: floor(wantedK * PICK_STRONGEST_RATIO))];
end

function [X, Y] = preprocess(mergedPeaksMin)
    X = mergedPeaksMin.seqs;
    L = size(X, 2);
    Y = (mergedPeaksMin.overlaps > 0) * (1:size(mergedPeaksMin.overlaps, 2))';
    [overlaps, seqInd] = sortrows(mergedPeaksMin.overlaps);
    X = X(seqInd, round(2 * L / 5):round(3 * L / 5));
    Y = Y(seqInd);
end
